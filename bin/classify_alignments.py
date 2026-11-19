#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Classify ComRepka-like flank alignments from minimap2 PAF.

Expected PAF command example:
  minimap2 -x asm10 -t 64 -c --cs=long --secondary=yes -N 20 -p 0.8 REF.mmi flanks.fa > flanks.paf

Expected query header format from the ComRepka-like flank builder:
  >group_id|chrom|left|L=...|left:start-end|block:start-end|...|overhang:50
  >group_id|chrom|right|L=...|right:start-end|block:start-end|...|overhang:50

Expected flank table columns:
  group_id chrom group_start group_end blocked_start blocked_end
  left_flank_start left_flank_end right_flank_start right_flank_end ...

Main output: one row per reference/main-genome locus with a top-level status
and detailed diagnostic fields.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Iterable
from tqdm import tqdm


# -----------------------------
# Data classes
# -----------------------------
@dataclass
class FlankInfo:
    group_id: str
    chrom: str
    group_start: int
    group_end: int
    left_flank_start: int
    left_flank_end: int
    right_flank_start: int
    right_flank_end: int
    repeat_ids_inside: str = ""

    @property
    def repeat_len(self) -> int:
        return self.group_end - self.group_start

    @property
    def left_overhang_len(self) -> int:
        # left query ends inside repeat: left_flank_end = group_start + overhang
        return max(0, self.left_flank_end - self.group_start)

    @property
    def right_overhang_len(self) -> int:
        # right query starts inside repeat: right_flank_start = group_end - overhang
        return max(0, self.group_end - self.right_flank_start)


@dataclass
class PAFAln:
    qname: str
    group_id: str
    side: str  # left/right
    qlen: int
    qstart: int
    qend: int
    strand: str
    tname: str
    tlen: int
    tstart: int
    tend: int
    matches: int
    block_len: int
    mapq: int
    tags: Dict[str, str] = field(default_factory=dict)

    @property
    def identity(self) -> float:
        if self.block_len <= 0:
            return 0.0
        return self.matches / self.block_len

    @property
    def qspan(self) -> int:
        return max(0, self.qend - self.qstart)

    @property
    def qcov(self) -> float:
        if self.qlen <= 0:
            return 0.0
        return self.qspan / self.qlen

    @property
    def tp(self) -> str:
        return self.tags.get("tp", "")

    @property
    def divergence(self) -> str:
        return self.tags.get("de", "")

    @property
    def aln_rank_score(self) -> float:
        # Used only for pruning stored alignments per side.
        primary_bonus = 0.05 if self.tp == "P" else 0.0
        mapq_bonus = min(self.mapq, 60) / 600.0
        return (0.55 * self.qcov) + (0.45 * self.identity) + primary_bonus + mapq_bonus


@dataclass
class GroupAlignments:
    left: List[PAFAln] = field(default_factory=list)
    right: List[PAFAln] = field(default_factory=list)
    left_raw: int = 0
    right_raw: int = 0
    left_filtered: int = 0
    right_filtered: int = 0


@dataclass
class PairEval:
    left: PAFAln
    right: PAFAln
    delta: int
    pair_span: int
    observed_interval_len: int
    length_delta: int
    tsd_len: Optional[int]
    tsd_start: Optional[int]
    tsd_end: Optional[int]
    left_overhang_supported: bool
    right_overhang_supported: bool
    candidate_status: str
    status_reason: str
    score: float


# -----------------------------
# Parsers
# -----------------------------
def parse_int(value: str, field_name: str, line_no: int) -> int:
    try:
        return int(value)
    except ValueError:
        raise ValueError(f"Line {line_no}: cannot parse integer field {field_name}={value!r}")


def load_flank_table(path: str) -> "OrderedDict[str, FlankInfo]":
    out: "OrderedDict[str, FlankInfo]" = OrderedDict()
    with open(path, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = [
            "group_id", "chrom", "group_start", "group_end",
            "left_flank_start", "left_flank_end",
            "right_flank_start", "right_flank_end",
        ]
        missing = [c for c in required if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f"Missing required columns in flank table: {', '.join(missing)}")
        for i, row in enumerate(tqdm(reader, desc="Loading flank table", unit="rows"), start=2):
            gid = row["group_id"]
            if not gid:
                continue
            if gid in out:
                raise ValueError(f"Duplicate group_id in flank table at line {i}: {gid}")
            out[gid] = FlankInfo(
                group_id=gid,
                chrom=row["chrom"],
                group_start=parse_int(row["group_start"], "group_start", i),
                group_end=parse_int(row["group_end"], "group_end", i),
                left_flank_start=parse_int(row["left_flank_start"], "left_flank_start", i),
                left_flank_end=parse_int(row["left_flank_end"], "left_flank_end", i),
                right_flank_start=parse_int(row["right_flank_start"], "right_flank_start", i),
                right_flank_end=parse_int(row["right_flank_end"], "right_flank_end", i),
                repeat_ids_inside=row.get("repeat_ids_inside", ""),
            )
    return out


def parse_query_name(qname: str) -> Tuple[str, Optional[str], Optional[int]]:
    """
    Returns (group_id, side, header_overhang_len).
    side is parsed from a literal token 'left'/'right'.
    """
    parts = qname.split("|")
    group_id = parts[0]
    side = None
    header_ov = None
    for p in parts[1:]:
        if p == "left" or p == "right":
            side = p
        elif p.startswith("overhang:"):
            try:
                header_ov = int(p.split(":", 1)[1])
            except ValueError:
                header_ov = None
    return group_id, side, header_ov


def parse_paf_tags(fields: List[str]) -> Dict[str, str]:
    tags: Dict[str, str] = {}
    for f in fields:
        # tag:type:value ; value may contain ':' so split max 2
        parts = f.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def parse_paf_line(line: str, line_no: int) -> Optional[PAFAln]:
    if not line.strip() or line.startswith("#"):
        return None
    fields = line.rstrip("\n").split()
    if len(fields) < 12:
        raise ValueError(f"Line {line_no}: PAF has <12 columns")

    qname = fields[0]
    group_id, side, _ = parse_query_name(qname)
    if side not in {"left", "right"}:
        # Unknown query naming convention; skip.
        return None

    return PAFAln(
        qname=qname,
        group_id=group_id,
        side=side,
        qlen=parse_int(fields[1], "qlen", line_no),
        qstart=parse_int(fields[2], "qstart", line_no),
        qend=parse_int(fields[3], "qend", line_no),
        strand=fields[4],
        tname=fields[5],
        tlen=parse_int(fields[6], "tlen", line_no),
        tstart=parse_int(fields[7], "tstart", line_no),
        tend=parse_int(fields[8], "tend", line_no),
        matches=parse_int(fields[9], "matches", line_no),
        block_len=parse_int(fields[10], "block_len", line_no),
        mapq=parse_int(fields[11], "mapq", line_no),
        tags=parse_paf_tags(fields[12:]),
    )


# -----------------------------
# Filtering and support checks
# -----------------------------
def passes_alignment_filters(aln: PAFAln, args: argparse.Namespace) -> bool:
    if aln.identity < args.min_identity:
        return False
    if aln.qcov < args.min_query_cover:
        return False
    if aln.mapq < args.min_mapq:
        return False
    if args.primary_only and aln.tp and aln.tp != "P":
        return False
    return True


def add_alignment(group: GroupAlignments, aln: PAFAln, args: argparse.Namespace, filtered: bool) -> None:
    if aln.side == "left":
        group.left_raw += 1
        if filtered:
            group.left_filtered += 1
            group.left.append(aln)
            prune_alignments(group.left, args.max_hits_per_side)
    elif aln.side == "right":
        group.right_raw += 1
        if filtered:
            group.right_filtered += 1
            group.right.append(aln)
            prune_alignments(group.right, args.max_hits_per_side)


def prune_alignments(items: List[PAFAln], max_n: int) -> None:
    if max_n <= 0:
        return
    if len(items) > max_n * 2:
        items.sort(key=lambda x: x.aln_rank_score, reverse=True)
        del items[max_n:]


def final_prune(items: List[PAFAln], max_n: int) -> None:
    if max_n > 0 and len(items) > max_n:
        items.sort(key=lambda x: x.aln_rank_score, reverse=True)
        del items[max_n:]


def overhang_supported(aln: PAFAln, side: str, ov_len: int, args: argparse.Namespace) -> bool:
    """
    Query-coordinate support for the repeat overhang.

    left query:  [flank][overhang]  => overhang interval [qlen-ov, qlen)
    right query: [overhang][flank]  => overhang interval [0, ov)

    This checks that the aligned query span covers enough of that interval and
    reaches the corresponding query edge within --overhang-edge-slop.
    """
    if ov_len <= 0:
        return False
    min_cov = min(args.min_overhang_aligned, ov_len)

    if side == "left":
        ov_s, ov_e = max(0, aln.qlen - ov_len), aln.qlen
        covered = max(0, min(aln.qend, ov_e) - max(aln.qstart, ov_s))
        edge_ok = aln.qend >= aln.qlen - args.overhang_edge_slop
        return covered >= min_cov and edge_ok

    if side == "right":
        ov_s, ov_e = 0, min(ov_len, aln.qlen)
        covered = max(0, min(aln.qend, ov_e) - max(aln.qstart, ov_s))
        edge_ok = aln.qstart <= args.overhang_edge_slop
        return covered >= min_cov and edge_ok

    return False


# -----------------------------
# Pair geometry and classification
# -----------------------------
def eval_pair(left: PAFAln, right: PAFAln, info: FlankInfo, args: argparse.Namespace) -> Optional[PairEval]:
    if left.tname != right.tname:
        return None
    if left.strand != right.strand:
        return None
    if left.strand not in {"+", "-"}:
        return None

    # Correct order and delta between repeat-proximal target edges.
    if left.strand == "+":
        # left interval should start before right interval. Overlap near the boundary is allowed.
        if left.tstart > right.tstart:
            return None
        delta = right.tstart - left.tend
        pair_span = max(left.tend, right.tend) - min(left.tstart, right.tstart)
        if delta <= 0:
            tsd_len = -delta
            tsd_start = right.tstart
            tsd_end = left.tend
        else:
            tsd_len = None
            tsd_start = None
            tsd_end = None
    else:
        # On '-' strand, the right flank maps before the left flank in target coordinates.
        if right.tstart > left.tstart:
            return None
        delta = left.tstart - right.tend
        pair_span = max(left.tend, right.tend) - min(left.tstart, right.tstart)
        if delta <= 0:
            tsd_len = -delta
            tsd_start = left.tstart
            tsd_end = right.tend
        else:
            tsd_len = None
            tsd_start = None
            tsd_end = None

    left_ov = info.left_overhang_len
    right_ov = info.right_overhang_len
    left_oh = overhang_supported(left, "left", left_ov, args)
    right_oh = overhang_supported(right, "right", right_ov, args)

    observed_interval_len = delta + left_ov + right_ov
    length_delta = observed_interval_len - info.repeat_len

    # Do not treat very distant pairs as local orthologous candidates.
    pair_dist_limit = max(args.max_pair_distance, info.repeat_len + args.length_tolerance_bp)
    if abs(delta) > pair_dist_limit:
        return None

    tol = max(args.length_tolerance_bp, int(math.ceil(args.length_tolerance_frac * max(1, info.repeat_len))))

    if left_oh and right_oh:
        if abs(length_delta) <= tol:
            cand_status = "occupied_expected_length"
            reason = "both_repeat_overhangs_supported_and_observed_length_within_tolerance"
            status_bonus = 35.0
        else:
            cand_status = "occupied_different_length"
            reason = "both_repeat_overhangs_supported_but_observed_length_outside_tolerance"
            status_bonus = 25.0
    elif (not left_oh) and (not right_oh):
        if delta <= 0 and tsd_len is not None and 0 <= tsd_len <= args.max_tsd:
            cand_status = "empty_with_TSD"
            reason = "both_repeat_overhangs_not_supported_and_flank_alignments_share_short_target_interval"
            status_bonus = 35.0
        elif delta > 0 and delta <= args.max_empty_gap:
            cand_status = "empty_without_TSD"
            reason = "both_repeat_overhangs_not_supported_and_flanks_have_small_gap_without_TSD_overlap"
            status_bonus = 20.0
        else:
            cand_status = "local_indel_or_rearranged"
            reason = "both_repeat_overhangs_not_supported_but_flank_geometry_is_not_clean_empty_locus"
            status_bonus = 5.0
    else:
        cand_status = "local_indel_or_rearranged"
        reason = "only_one_repeat_overhang_supported"
        status_bonus = 5.0

    avg_identity = (left.identity + right.identity) / 2.0
    avg_qcov = (left.qcov + right.qcov) / 2.0
    primary_bonus = 0.0
    if left.tp == "P":
        primary_bonus += 3.0
    if right.tp == "P":
        primary_bonus += 3.0

    # Penalize very large unexpected separation, but keep the penalty modest.
    distance_penalty = min(abs(delta) / max(1.0, args.max_pair_distance) * 10.0, 10.0)

    score = (45.0 * avg_identity) + (45.0 * avg_qcov) + status_bonus + primary_bonus - distance_penalty

    return PairEval(
        left=left,
        right=right,
        delta=delta,
        pair_span=pair_span,
        observed_interval_len=observed_interval_len,
        length_delta=length_delta,
        tsd_len=tsd_len,
        tsd_start=tsd_start,
        tsd_end=tsd_end,
        left_overhang_supported=left_oh,
        right_overhang_supported=right_oh,
        candidate_status=cand_status,
        status_reason=reason,
        score=score,
    )


def classify_group(info: FlankInfo, galn: GroupAlignments, args: argparse.Namespace) -> Dict[str, str]:
    final_prune(galn.left, args.max_hits_per_side)
    final_prune(galn.right, args.max_hits_per_side)

    base = {
        "group_id": info.group_id,
        "ref_chrom": info.chrom,
        "ref_start": str(info.group_start),
        "ref_end": str(info.group_end),
        "repeat_len": str(info.repeat_len),
        "repeat_ids_inside": info.repeat_ids_inside,
        "left_raw_hits": str(galn.left_raw),
        "right_raw_hits": str(galn.right_raw),
        "left_filtered_hits": str(galn.left_filtered),
        "right_filtered_hits": str(galn.right_filtered),
    }

    if galn.left_filtered == 0 and galn.right_filtered == 0:
        return finish_no_pair(base, "no_flank_hits", "no_filtered_left_or_right_alignments", args)
    if galn.left_filtered == 0:
        return finish_no_pair(base, "no_left_hit", "no_filtered_left_alignments", args)
    if galn.right_filtered == 0:
        return finish_no_pair(base, "no_right_hit", "no_filtered_right_alignments", args)

    pairs: List[PairEval] = []
    same_target_order_ignored_distance = 0

    for l in galn.left:
        for r in galn.right:
            # Count rough same-locus attempts even if distance filtering later rejects.
            if l.tname == r.tname and l.strand == r.strand:
                if (l.strand == "+" and l.tstart <= r.tstart) or (l.strand == "-" and r.tstart <= l.tstart):
                    same_target_order_ignored_distance += 1
            pe = eval_pair(l, r, info, args)
            if pe is not None:
                pairs.append(pe)

    if not pairs:
        if same_target_order_ignored_distance > 0:
            return finish_no_pair(base, "discordant_pair", "same_contig_strand_order_exists_but_distance_exceeds_limit", args)
        return finish_no_pair(base, "no_valid_lr_pair", "no_same_contig_strand_ordered_left_right_pair", args)

    pairs.sort(key=lambda p: p.score, reverse=True)
    best = pairs[0]
    second = pairs[1] if len(pairs) > 1 else None

    status = best.candidate_status
    status_reason = best.status_reason
    if second is not None and (best.score - second.score) < args.ambiguous_score_delta:
        status = "ambiguous_lr_pair"
        status_reason = (
            f"best_pair_score_close_to_second_pair:"
            f"{best.score:.3f}-{second.score:.3f}<delta_{args.ambiguous_score_delta}"
        )

    row = dict(base)
    row.update({
        "status": status,
        "best_candidate_status": best.candidate_status,
        "status_reason": status_reason,
        "candidate_pair_count": str(len(pairs)),
        "best_pair_score": f"{best.score:.3f}",
        "second_pair_score": "" if second is None else f"{second.score:.3f}",
        "target_contig": best.left.tname,
        "target_strand": best.left.strand,
        "left_target_start": str(best.left.tstart),
        "left_target_end": str(best.left.tend),
        "right_target_start": str(best.right.tstart),
        "right_target_end": str(best.right.tend),
        "left_query_start": str(best.left.qstart),
        "left_query_end": str(best.left.qend),
        "right_query_start": str(best.right.qstart),
        "right_query_end": str(best.right.qend),
        "left_query_len": str(best.left.qlen),
        "right_query_len": str(best.right.qlen),
        "left_identity": f"{best.left.identity:.6f}",
        "right_identity": f"{best.right.identity:.6f}",
        "left_query_cover": f"{best.left.qcov:.6f}",
        "right_query_cover": f"{best.right.qcov:.6f}",
        "left_mapq": str(best.left.mapq),
        "right_mapq": str(best.right.mapq),
        "left_tp": best.left.tp,
        "right_tp": best.right.tp,
        "left_divergence": best.left.divergence,
        "right_divergence": best.right.divergence,
        "left_overhang_len": str(info.left_overhang_len),
        "right_overhang_len": str(info.right_overhang_len),
        "left_overhang_supported": str(best.left_overhang_supported),
        "right_overhang_supported": str(best.right_overhang_supported),
        "delta": str(best.delta),
        "observed_interval_len": str(best.observed_interval_len),
        "length_delta": str(best.length_delta),
        "pair_span": str(best.pair_span),
        "tsd_len": "" if best.tsd_len is None else str(best.tsd_len),
        "tsd_target_start": "" if best.tsd_start is None else str(best.tsd_start),
        "tsd_target_end": "" if best.tsd_end is None else str(best.tsd_end),
        "left_qname": best.left.qname,
        "right_qname": best.right.qname,
    })

    confidence, low_conf, reasons = assess_confidence(row, best, status, args)
    row["confidence"] = confidence
    row["low_confidence"] = str(low_conf)
    row["confidence_reasons"] = ";".join(reasons)
    return fill_missing_output_columns(row)


def finish_no_pair(base: Dict[str, str], status: str, reason: str, args: argparse.Namespace) -> Dict[str, str]:
    row = dict(base)
    row.update({
        "status": status,
        "best_candidate_status": "",
        "status_reason": reason,
        "candidate_pair_count": "0",
        "best_pair_score": "",
        "second_pair_score": "",
        "confidence": "low",
        "low_confidence": "True",
        "confidence_reasons": status,
    })
    return fill_missing_output_columns(row)


def assess_confidence(row: Dict[str, str], best: PairEval, status: str, args: argparse.Namespace) -> Tuple[str, bool, List[str]]:
    severe: List[str] = []
    medium: List[str] = []

    if status in {"no_flank_hits", "no_left_hit", "no_right_hit", "no_valid_lr_pair", "discordant_pair", "ambiguous_lr_pair"}:
        severe.append(status)

    if best.left.identity < args.low_identity or best.right.identity < args.low_identity:
        severe.append("identity_below_low_identity")
    elif best.left.identity < args.high_identity or best.right.identity < args.high_identity:
        medium.append("identity_below_high_identity")

    if best.left.qcov < args.low_query_cover or best.right.qcov < args.low_query_cover:
        severe.append("query_cover_below_low_query_cover")
    elif best.left.qcov < args.high_query_cover or best.right.qcov < args.high_query_cover:
        medium.append("query_cover_below_high_query_cover")

    if args.high_mapq >= 0 and (best.left.mapq < args.high_mapq or best.right.mapq < args.high_mapq):
        medium.append("mapq_below_high_mapq")

    if status in {"empty_without_TSD", "occupied_different_length", "local_indel_or_rearranged"}:
        medium.append(f"status_{status}")

    if severe:
        return "low", True, severe + medium
    if medium:
        return "medium", False, medium
    return "high", False, []


OUTPUT_COLUMNS = [
    "group_id", "ref_chrom", "ref_start", "ref_end", "repeat_len", "repeat_ids_inside",
    "status", "best_candidate_status", "confidence", "low_confidence", "confidence_reasons", "status_reason",
    "left_raw_hits", "right_raw_hits", "left_filtered_hits", "right_filtered_hits", "candidate_pair_count",
    "best_pair_score", "second_pair_score",
    "target_contig", "target_strand",
    "left_target_start", "left_target_end", "right_target_start", "right_target_end",
    "left_query_start", "left_query_end", "right_query_start", "right_query_end",
    "left_query_len", "right_query_len",
    "left_identity", "right_identity", "left_query_cover", "right_query_cover",
    "left_mapq", "right_mapq", "left_tp", "right_tp", "left_divergence", "right_divergence",
    "left_overhang_len", "right_overhang_len", "left_overhang_supported", "right_overhang_supported",
    "delta", "observed_interval_len", "length_delta", "pair_span",
    "tsd_len", "tsd_target_start", "tsd_target_end",
    "left_qname", "right_qname",
]


def fill_missing_output_columns(row: Dict[str, str]) -> Dict[str, str]:
    for c in OUTPUT_COLUMNS:
        row.setdefault(c, "")
    return row


# -----------------------------
# I/O modes
# -----------------------------
def run_stream_mode(flank_infos: "OrderedDict[str, FlankInfo]", args: argparse.Namespace) -> None:
    emitted = set()
    closed_groups = set()

    with open(args.out_table, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()

        current_gid: Optional[str] = None
        current = GroupAlignments()

        def flush_current() -> None:
            nonlocal current_gid, current
            if current_gid is None:
                return
            info = flank_infos.get(current_gid)
            if info is not None:
                writer.writerow(classify_group(info, current, args))
                emitted.add(current_gid)
            elif not args.skip_unknown_groups:
                # Unknown group: cannot compute repeat length/overhang from table, so skip classification safely.
                sys.stderr.write(f"[WARN] group_id not found in flank table, skipped: {current_gid}\n")
            closed_groups.add(current_gid)
            current_gid = None
            current = GroupAlignments()

        with open(args.paf, "r") as paf_fh:
            for line_no, line in enumerate(tqdm(paf_fh, desc="Reading PAF", unit="lines"), start=1):
                aln = parse_paf_line(line, line_no)
                if aln is None:
                    continue
                if current_gid is None:
                    current_gid = aln.group_id
                elif aln.group_id != current_gid:
                    flush_current()
                    if aln.group_id in closed_groups:
                        raise ValueError(f"Non-contiguous PAF group {aln.group_id}; rerun with --mode collect")
                    current_gid = aln.group_id
                    current = GroupAlignments()

                filtered = passes_alignment_filters(aln, args)
                add_alignment(current, aln, args, filtered)

        flush_current()

        # Emit groups absent from PAF entirely.
        for gid, info in tqdm(flank_infos.items(), desc="Writing absent PAF groups", unit="groups"):
            if gid not in emitted:
                writer.writerow(classify_group(info, GroupAlignments(), args))


def run_collect_mode(flank_infos: "OrderedDict[str, FlankInfo]", args: argparse.Namespace) -> None:
    groups: Dict[str, GroupAlignments] = {}

    with open(args.paf, "r") as paf_fh:
        for line_no, line in enumerate(tqdm(paf_fh, desc="Reading PAF", unit="lines"), start=1):
            aln = parse_paf_line(line, line_no)
            if aln is None:
                continue
            if aln.group_id not in flank_infos:
                if not args.skip_unknown_groups:
                    sys.stderr.write(f"[WARN] group_id not found in flank table, skipped: {aln.group_id}\n")
                continue
            galn = groups.setdefault(aln.group_id, GroupAlignments())
            filtered = passes_alignment_filters(aln, args)
            add_alignment(galn, aln, args, filtered)

    with open(args.out_table, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gid, info in tqdm(flank_infos.items(), desc="Classifying groups", unit="groups"):
            writer.writerow(classify_group(info, groups.get(gid, GroupAlignments()), args))


# -----------------------------
# CLI
# -----------------------------
def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Classify ComRepka-like minimap2 PAF flank alignments into occupied/empty/rearranged statuses."
    )
    ap.add_argument("--paf", required=True, help="minimap2 PAF file for ComRepka-like left/right flank queries")
    ap.add_argument("--flank-table", required=True, help="TSV table produced by the flank-preparation script")
    ap.add_argument("--out-table", required=True, help="Output TSV classification table")

    ap.add_argument("--mode", choices=["stream", "collect"], default="stream",
                    help="stream assumes PAF is grouped by query/group, as minimap2 normally outputs; collect is safer but uses more RAM")
    ap.add_argument("--skip-unknown-groups", action="store_true", default=True,
                    help="Skip PAF records whose group_id is absent from flank table")

    # Core thresholds
    ap.add_argument("--max-tsd", type=int, default=50,
                    help="Maximum target site duplication length allowed for empty_with_TSD")
    ap.add_argument("--max-empty-gap", type=int, default=50,
                    help="Maximum small positive gap for empty_without_TSD")
    ap.add_argument("--length-tolerance-bp", type=int, default=100,
                    help="Absolute length tolerance for occupied_expected_length")
    ap.add_argument("--length-tolerance-frac", type=float, default=0.10,
                    help="Fractional length tolerance for occupied_expected_length; max(bp, frac*repeat_len) is used")
    ap.add_argument("--max-pair-distance", type=int, default=10000,
                    help="Maximum absolute distance between repeat-proximal flank edges for a local pair; automatically raised to repeat_len + tolerance_bp if needed")

    # Alignment filters
    ap.add_argument("--min-identity", type=float, default=0.80,
                    help="Minimum PAF identity = matches / alignment_block_length")
    ap.add_argument("--min-query-cover", type=float, default=0.20,
                    help="Minimum query span coverage = (qend-qstart)/qlen")
    ap.add_argument("--min-mapq", type=int, default=0,
                    help="Minimum PAF MAPQ")
    ap.add_argument("--primary-only", action="store_true",
                    help="Use only alignments tagged tp:A:P")
    ap.add_argument("--max-hits-per-side", type=int, default=20,
                    help="Keep only top N filtered alignments per side per group before pairing")

    # Overhang support
    ap.add_argument("--min-overhang-aligned", type=int, default=45,
                    help="Minimum query bases from the repeat-overhang interval that must be covered")
    ap.add_argument("--overhang-edge-slop", type=int, default=5,
                    help="Allowed query-edge clipping while still calling overhang support")

    # Ambiguity and confidence
    ap.add_argument("--ambiguous-score-delta", type=float, default=5.0,
                    help="If best_pair_score - second_pair_score is below this value, status becomes ambiguous_lr_pair")
    ap.add_argument("--high-identity", type=float, default=0.90,
                    help="Identity threshold for high confidence")
    ap.add_argument("--low-identity", type=float, default=0.85,
                    help="Identity threshold below which confidence is low")
    ap.add_argument("--high-query-cover", type=float, default=0.50,
                    help="Query coverage threshold for high confidence")
    ap.add_argument("--low-query-cover", type=float, default=0.30,
                    help="Query coverage threshold below which confidence is low")
    ap.add_argument("--high-mapq", type=int, default=-1,
                    help="Optional MAPQ threshold for high confidence; -1 disables MAPQ confidence check")

    return ap


def main() -> None:
    args = build_argparser().parse_args()

    flank_infos = load_flank_table(args.flank_table)
    if not flank_infos:
        raise SystemExit("No loci loaded from flank table")

    sys.stderr.write(f"[INFO] Loaded loci from flank table: {len(flank_infos)}\n")
    sys.stderr.write(f"[INFO] Mode: {args.mode}\n")

    if args.mode == "stream":
        run_stream_mode(flank_infos, args)
    else:
        run_collect_mode(flank_infos, args)

    sys.stderr.write(f"[OK] Wrote: {args.out_table}\n")


if __name__ == "__main__":
    main()
