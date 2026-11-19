#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HAL anchor-based validation script (single-file).

High-level idea:
  1) For each group_id, build anchor windows inside left/right flanks in SOURCE genome.
  2) Lift anchors to TARGET genome using halLiftover (batched).
  3) For each group_id, form candidate gap intervals G from pairs:
        (left-anchor hit) + (right-anchor hit)
     Candidate validity (STRICT orientation, Variant A):
        - same chromosome
        - order consistent with pred_strand
        - g_end > g_start
  4) Compare each candidate G to predicted interval P from minimap2:
        overlap_bp = |G ∩ P|
        pred_recall = overlap_bp / |P|
        IoU = |G ∩ P| / |G ∪ P|   (IoU=0 if different chromosomes)
        len_ratio = max(|G|, orig_len) / max(1, min(|G|, orig_len))
           where orig_len = right_flank_start - left_flank_end (SOURCE inner-boundary gap length)
  5) Status logic:
        - HAL_UNSUPPORTED:
            no valid candidates at all
        - HAL_DISAGREE:
            valid candidates exist but Support is empty (max pred_recall < R)
        - HAL_CONFIRMED_UNIQUE:
            Support non-empty AND no "real alternative other area" candidates
        - HAL_AMBIGUOUS:
            Support non-empty AND exists at least one "real alternative other area"

     Definitions:
        Support:   pred_recall >= R
        Alt:       pred_recall < R
        AltOther:  Alt and IoU < I
        AltOtherReal: AltOther and len_ratio <= 2.0   (length comparability to orig_len)

Outputs:
  - out_main: one row per group_id with status + best SUPPORT candidate (if exists)
  - out_amb : candidates for AMBIGUOUS groups (all Support + all AltOtherReal)
  - out_dis : candidates for DISAGREE groups (all valid candidates with pred_recall < R)
  - out_other: optional debug (all valid candidates on chromosomes != pred_chr_hal)
  - out_unsup: debug for UNSUPPORTED (per offset-combo, with hit counts)
"""

import argparse
import csv
import json
import os
import sys
import tempfile
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

# -------- CSV safety --------
try:
    csv.field_size_limit(sys.maxsize)
except Exception:
    pass

# -------- Optional progress bar --------
try:
    from tqdm import tqdm
except Exception:
    tqdm = None


@dataclass(frozen=True)
class Interval:
    chrom: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return max(0, self.end - self.start)


def run_cmd(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if check and p.returncode != 0:
        raise RuntimeError(
            "Command failed\n"
            f"cmd: {' '.join(cmd)}\n"
            f"returncode: {p.returncode}\n"
            f"stdout:\n{p.stdout}\n"
            f"stderr:\n{p.stderr}\n"
        )
    return p


def load_json_map(path: Optional[str]) -> Dict[str, str]:
    if not path:
        return {}
    with open(path, "r") as f:
        return json.load(f)


def map_chrom(name: str, mapping: Dict[str, str]) -> str:
    return mapping.get(name, name)


def overlap_len(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def iou_1d(a0: int, a1: int, b0: int, b1: int) -> float:
    """IoU for 1D half-open intervals [a0,a1), [b0,b1)."""
    ov = overlap_len(a0, a1, b0, b1)
    la = max(0, a1 - a0)
    lb = max(0, b1 - b0)
    union = la + lb - ov
    if union <= 0:
        return 0.0
    return ov / union


def clamp_anchor_within_flank(
    flank_start: int,
    flank_end: int,
    anchor_start: int,
    anchor_end: int,
) -> Optional[Tuple[int, int]]:
    # Anchor must be fully inside flank
    if anchor_start < flank_start or anchor_end > flank_end:
        return None
    if anchor_end <= anchor_start:
        return None
    return (anchor_start, anchor_end)


def parse_offsets(s: str) -> List[int]:
    out = []
    for x in s.split(","):
        x = x.strip()
        if not x:
            continue
        out.append(int(x))
    return out


def load_flank_table(flank_table_tsv: str) -> Dict[str, Dict[str, int]]:
    """
    Expected columns:
      group_id, chrom, left_flank_start,left_flank_end,right_flank_start,right_flank_end
    """
    d: Dict[str, Dict[str, int]] = {}
    with open(flank_table_tsv, "r") as f:
        r = csv.DictReader(f, delimiter="\t")
        need = [
            "group_id",
            "chrom",
            "left_flank_start",
            "left_flank_end",
            "right_flank_start",
            "right_flank_end",
        ]
        for c in need:
            if c not in (r.fieldnames or []):
                raise RuntimeError(f"flank-table missing column: {c}")
        for row in r:
            gid = row["group_id"]
            d[gid] = {
                "src_chrom": row["chrom"],
                "left_flank_start": int(row["left_flank_start"]),
                "left_flank_end": int(row["left_flank_end"]),
                "right_flank_start": int(row["right_flank_start"]),
                "right_flank_end": int(row["right_flank_end"]),
            }
    return d


def load_pairs(pairs_tsv: str) -> Dict[str, Dict[str, object]]:
    """
    Load predicted gap intervals P on the TARGET genome.

    Supports two input formats (auto-detected by column names):

      (A) Legacy format:
          group_id, target_chr, target_strand, target_gap_start, target_gap_end

      (B) New ComRepka-style table (wide table):
          group_id, target_contig, target_strand, left_target_end, right_target_start

        Here the predicted gap is defined by the inner boundaries between flanks in target:
          pred_start = min(left_target_end, right_target_start)
          pred_end   = max(left_target_end, right_target_start)

          + strand: left_target_end -> right_target_start
          - strand: right_target_end -> left_target_start

    Notes:
      - Coordinates are assumed 0-based, half-open [start,end).
      - If required prediction fields are missing/empty for a row, the row is still
        included but pred_start==pred_end==0 so that downstream logic can emit
        HAL_UNSUPPORTED with a diagnostic note.
    """
    out: Dict[str, Dict[str, object]] = {}
    with open(pairs_tsv, "r") as f:
        r = csv.DictReader(f, delimiter="	")
        fields = set(r.fieldnames or [])

        # Decide format by presence of key columns
        has_legacy = {"target_chr", "target_gap_start", "target_gap_end"}.issubset(fields)
        has_new = {
                        "target_contig",
                        "left_target_start",
                        "left_target_end",
                        "right_target_start",
                        "right_target_end",
                    }.issubset(fields)

        if not has_legacy and not has_new:
            raise RuntimeError(
                "pairs-tsv has unrecognized columns. Expected either legacy columns "
                "(target_chr/target_gap_start/target_gap_end) or new columns "
                "(target_contig/left_target_end/right_target_start)."
            )

        # target_strand and group_id are required in both formats
        for req in ["group_id", "target_strand"]:
            if req not in fields:
                raise RuntimeError(f"pairs-tsv missing column: {req}")

        for row in r:
            gid = row.get("group_id", "")
            if not gid:
                continue

            pred_strand = str(row.get("target_strand", "")).strip()
            if not pred_strand:
                pred_strand = "+"  # fallback, but keep row

            if has_legacy:
                pred_chr = str(row.get("target_chr", "")).strip()
                s_raw = str(row.get("target_gap_start", "")).strip()
                e_raw = str(row.get("target_gap_end", "")).strip()

                if pred_chr and s_raw and e_raw:
                    s = int(s_raw)
                    e = int(e_raw)
                    pred_start, pred_end = (s, e) if s <= e else (e, s)
                else:
                    # Missing prediction -> keep row but with empty coords
                    pred_chr = pred_chr or ""
                    pred_start = 0
                    pred_end = 0

            else:
                pred_chr = str(row.get("target_contig", "")).strip()
                
                if pred_strand == "+":
                    a_raw = str(row.get("left_target_end", "")).strip()
                    b_raw = str(row.get("right_target_start", "")).strip()
                else:
                    a_raw = str(row.get("right_target_end", "")).strip()
                    b_raw = str(row.get("left_target_start", "")).strip()
                
                if pred_chr and a_raw and b_raw:
                    a = int(a_raw)
                    b = int(b_raw)
                    pred_start, pred_end = (a, b) if a <= b else (b, a)
                else:
                    pred_chr = pred_chr or ""
                    pred_start = 0
                    pred_end = 0

            out[gid] = {
                "pred_chr_in": pred_chr,
                "pred_strand": pred_strand,
                "pred_start": int(pred_start),
                "pred_end": int(pred_end),
            }

    return out


def write_bed4(path: str, rows: List[Tuple[str, int, int, str]]) -> None:
    with open(path, "w") as w:
        for chrom, s, e, name in rows:
            w.write(f"{chrom}\t{s}\t{e}\t{name}\n")


def parse_bed4(path: str) -> List[Tuple[str, int, int, str]]:
    out = []
    with open(path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            t = line.rstrip("\n").split("\t")
            if len(t) < 4:
                continue
            out.append((t[0], int(t[1]), int(t[2]), t[3]))
    return out


def hal_liftover_batch(
    hal: str,
    src_genome: str,
    dst_genome: str,
    bed_in: str,
    bed_out: str,
    executable: str = "halLiftover",
) -> None:
    run_cmd([executable, hal, src_genome, bed_in, dst_genome, bed_out], check=True)


def build_anchor_requests(
    flank_meta: Dict[str, Dict[str, int]],
    pairs: Dict[str, Dict[str, object]],
    src_chrom_map: Dict[str, str],
    anchor_len: int,
    anchor_offsets: List[int],
) -> List[Tuple[str, int, int, str]]:
    """
    Returns BED4 rows: (chrom, start, end, name)
    name encodes: gid|side|offset
    """
    req: List[Tuple[str, int, int, str]] = []
    for gid in pairs.keys():
        if gid not in flank_meta:
            continue
        fm = flank_meta[gid]
        src_chrom_refseq = fm["src_chrom"]
        src_chrom_hal = map_chrom(src_chrom_refseq, src_chrom_map)

        lf0, lf1 = fm["left_flank_start"], fm["left_flank_end"]     # left flank
        rf0, rf1 = fm["right_flank_start"], fm["right_flank_end"]   # right flank

        for off in anchor_offsets:
            # Left anchor: ends at (left_flank_end - off)
            a_end = lf1 - off
            a_start = a_end - anchor_len
            c = clamp_anchor_within_flank(lf0, lf1, a_start, a_end)
            if c is not None:
                s, e = c
                name = f"{gid}|L|{off}"
                req.append((src_chrom_hal, s, e, name))

            # Right anchor: starts at (right_flank_start + off)
            b_start = rf0 + off
            b_end = b_start + anchor_len
            c2 = clamp_anchor_within_flank(rf0, rf1, b_start, b_end)
            if c2 is not None:
                s2, e2 = c2
                name2 = f"{gid}|R|{off}"
                req.append((src_chrom_hal, s2, e2, name2))

    return req


def group_hits_by_name(liftover_rows: List[Tuple[str, int, int, str]]) -> Dict[str, List[Interval]]:
    hits: Dict[str, List[Interval]] = defaultdict(list)
    for chrom, s, e, name in liftover_rows:
        hits[name].append(Interval(chrom, s, e))
    return hits


def choose_best_support(cands: List[Dict[str, object]]) -> Dict[str, object]:
    """
    Best SUPPORT candidate ranking:
      - higher pred_recall
      - higher overlap_bp
      - closer gap_len to pred_len
    """
    return max(
        cands,
        key=lambda x: (
            float(x["pred_recall"]),
            int(x["overlap_bp"]),
            -abs(int(x["gap_len"]) - int(x["pred_len"])),
        ),
    )


def compute_len_ratio(gap_len: int, orig_len: int) -> float:
    """len_ratio = max(gap_len, orig_len) / max(1, min(gap_len, orig_len))."""
    a = int(max(0, gap_len))
    b = int(max(0, orig_len))
    mn = min(a, b)
    mx = max(a, b)
    return mx / max(1, mn)


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--flank-table", required=True)
    ap.add_argument("--pairs-tsv", required=True)

    ap.add_argument("--hal", required=True)
    ap.add_argument("--hal-liftover", default="halLiftover", help="Executable name or absolute path")
    ap.add_argument("--src-genome", required=True)
    ap.add_argument("--dst-genome", required=True)

    ap.add_argument("--src-chrom-map", default=None)
    ap.add_argument("--dst-chrom-map", default=None)

    ap.add_argument("--anchor-len", type=int, default=300)
    ap.add_argument("--anchor-offsets", default="0,400")

    # New status logic parameters
    ap.add_argument("--min-pred-recall", type=float, default=0.95)
    ap.add_argument("--alt-iou-threshold", type=float, default=0.2)
    ap.add_argument("--len-ratio-max", type=float, default=2.0)

    ap.add_argument("--batch-size", type=int, default=200000)
    ap.add_argument("--out-main", required=True)
    ap.add_argument("--out-amb", required=True)
    ap.add_argument("--out-dis", required=True)
    ap.add_argument("--out-other", required=True)
    ap.add_argument("--out-unsup", required=True)

    args = ap.parse_args()

    anchor_offsets = parse_offsets(args.anchor_offsets)
    if len(anchor_offsets) == 0:
        raise RuntimeError("anchor-offsets parsed empty")

    # All offset combinations (left_offset, right_offset)
    expected_combos = [(lo, ro) for lo in anchor_offsets for ro in anchor_offsets]

    src_map = load_json_map(args.src_chrom_map)
    dst_map = load_json_map(args.dst_chrom_map)

    flank_meta = load_flank_table(args.flank_table)
    pairs = load_pairs(args.pairs_tsv)

    # Build anchor BED requests (SOURCE coords, chrom mapped to HAL naming)
    req = build_anchor_requests(flank_meta, pairs, src_map, args.anchor_len, anchor_offsets)
    if len(req) == 0:
        raise RuntimeError("No anchor requests generated. Check flank-table/pairs overlap and coordinates.")

    # Liftover anchors in batches
    all_hits_by_name: Dict[str, List[Interval]] = defaultdict(list)

    with tempfile.TemporaryDirectory() as tmpdir:
        bed_in = os.path.join(tmpdir, "in.bed")
        bed_out = os.path.join(tmpdir, "out.bed")

        it = range(0, len(req), args.batch_size)
        if tqdm is not None:
            it = tqdm(it, desc="halLiftover batches", unit="batch")

        for i0 in it:
            chunk = req[i0:i0 + args.batch_size]
            write_bed4(bed_in, chunk)
            hal_liftover_batch(args.hal, args.src_genome, args.dst_genome, bed_in, bed_out, args.hal_liftover)

            rows = parse_bed4(bed_out)
            hits_chunk = group_hits_by_name(rows)
            for k, v in hits_chunk.items():
                all_hits_by_name[k].extend(v)

    # ---------------- Outputs ----------------

    main_fields = [
        "group_id",
        "pred_chr_in",
        "pred_chr_hal",
        "pred_strand",
        "pred_start",
        "pred_end",
        "pred_len",
        "orig_len",
        "hal_status",
        "n_valid_candidates",
        "n_support",
        "n_alt_other_real",
        "best_left_offset",
        "best_right_offset",
        "best_gap_chr",
        "best_gap_start",
        "best_gap_end",
        "best_gap_len",
        "best_overlap_bp",
        "best_pred_recall",
        "best_iou",
        "best_len_ratio",
        "R_min_pred_recall",
        "I_alt_iou_threshold",
        "len_ratio_max",
        "note",
    ]

    cand_fields = [
        "group_id",
        "pred_chr_hal",
        "pred_strand",
        "pred_start",
        "pred_end",
        "pred_len",
        "orig_len",
        "left_offset",
        "right_offset",
        "gap_chr",
        "gap_start",
        "gap_end",
        "gap_len",
        "overlap_bp",
        "pred_recall",
        "iou",
        "len_ratio",
        "class",
        "n_left_hits",
        "n_right_hits",
        "n_pairs_total",
    ]

    out_main = open(args.out_main, "w")
    out_amb = open(args.out_amb, "w")
    out_dis = open(args.out_dis, "w")
    out_other = open(args.out_other, "w")
    out_unsup = open(args.out_unsup, "w")

    wm = csv.DictWriter(out_main, delimiter="\t", fieldnames=main_fields)
    wa = csv.DictWriter(out_amb, delimiter="\t", fieldnames=cand_fields)
    wd = csv.DictWriter(out_dis, delimiter="\t", fieldnames=cand_fields)
    wo = csv.DictWriter(out_other, delimiter="\t", fieldnames=cand_fields)
    wu = csv.DictWriter(out_unsup, delimiter="\t", fieldnames=cand_fields)

    wm.writeheader()
    wa.writeheader()
    wd.writeheader()
    wo.writeheader()
    wu.writeheader()

    group_ids = list(pairs.keys())
    it2 = group_ids
    if tqdm is not None:
        it2 = tqdm(group_ids, desc="Validate groups", unit="group")

    for gid in it2:
        if gid not in flank_meta:
            # If flank metadata missing, we cannot derive anchors/orig_len reliably.
            # Treat it as unsupported (but still emit a minimal main record).
            wm.writerow({
                "group_id": gid,
                "pred_chr_in": str(pairs[gid]["pred_chr_in"]),
                "pred_chr_hal": map_chrom(str(pairs[gid]["pred_chr_in"]), dst_map),
                "pred_strand": str(pairs[gid]["pred_strand"]),
                "pred_start": int(pairs[gid]["pred_start"]),
                "pred_end": int(pairs[gid]["pred_end"]),
                "pred_len": max(0, int(pairs[gid]["pred_end"]) - int(pairs[gid]["pred_start"])),
                "orig_len": "",
                "hal_status": "HAL_UNSUPPORTED",
                "n_valid_candidates": 0,
                "n_support": 0,
                "n_alt_other_real": 0,
                "best_left_offset": "",
                "best_right_offset": "",
                "best_gap_chr": "",
                "best_gap_start": "",
                "best_gap_end": "",
                "best_gap_len": "",
                "best_overlap_bp": 0,
                "best_pred_recall": 0.0,
                "best_iou": 0.0,
                "best_len_ratio": "",
                "R_min_pred_recall": args.min_pred_recall,
                "I_alt_iou_threshold": args.alt_iou_threshold,
                "len_ratio_max": args.len_ratio_max,
                "note": "MISSING_FLANK_META",
            })
            continue

        # Predicted interval P from minimap2 (TARGET coords)
        pred_chr_in = str(pairs[gid]["pred_chr_in"])
        pred_chr_hal = map_chrom(pred_chr_in, dst_map)
        pred_strand = str(pairs[gid]["pred_strand"])
        pred_start = int(pairs[gid]["pred_start"])
        pred_end = int(pairs[gid]["pred_end"])
        pred_len = max(0, pred_end - pred_start)

        # orig_len (SOURCE inner-boundary gap length)
        fm = flank_meta[gid]
        orig_len = max(0, int(fm["right_flank_start"]) - int(fm["left_flank_end"]))

        # Collect ALL anchor hits (do NOT pre-filter by chromosome)
        left_hits_by_off: Dict[int, List[Interval]] = {}
        right_hits_by_off: Dict[int, List[Interval]] = {}

        for off in anchor_offsets:
            kL = f"{gid}|L|{off}"
            kR = f"{gid}|R|{off}"
            left_hits_by_off[off] = list(all_hits_by_name.get(kL, []))
            right_hits_by_off[off] = list(all_hits_by_name.get(kR, []))

        # If prediction is missing (e.g., upstream table has no target coordinates),
        # emit HAL_UNSUPPORTED and skip candidate evaluation.
        if (not pred_chr_in) or (pred_len == 0):
            wm.writerow({
                "group_id": gid,
                "pred_chr_in": pred_chr_in,
                "pred_chr_hal": pred_chr_hal,
                "pred_strand": pred_strand,
                "pred_start": pred_start,
                "pred_end": pred_end,
                "pred_len": pred_len,
                "orig_len": orig_len,
                "hal_status": "HAL_UNSUPPORTED",
                "n_valid_candidates": 0,
                "n_support": 0,
                "n_alt_other_real": 0,
                "best_left_offset": "",
                "best_right_offset": "",
                "best_gap_chr": "",
                "best_gap_start": "",
                "best_gap_end": "",
                "best_gap_len": "",
                "best_overlap_bp": 0,
                "best_pred_recall": 0.0,
                "best_iou": 0.0,
                "best_len_ratio": "",
                "R_min_pred_recall": args.min_pred_recall,
                "I_alt_iou_threshold": args.alt_iou_threshold,
                "len_ratio_max": args.len_ratio_max,
                "note": "MISSING_PREDICTION",
            })

            for lo, ro in expected_combos:
                Lhits = left_hits_by_off.get(lo, [])
                Rhits = right_hits_by_off.get(ro, [])
                n_pairs_total = len(Lhits) * len(Rhits)
                wu.writerow({
                    "group_id": gid,
                    "pred_chr_hal": pred_chr_hal,
                    "pred_strand": pred_strand,
                    "pred_start": pred_start,
                    "pred_end": pred_end,
                    "pred_len": pred_len,
                    "orig_len": orig_len,
                    "left_offset": lo,
                    "right_offset": ro,
                    "gap_chr": "",
                    "gap_start": "",
                    "gap_end": "",
                    "gap_len": "",
                    "overlap_bp": 0,
                    "pred_recall": 0.0,
                    "iou": 0.0,
                    "len_ratio": "",
                    "class": "NO_PREDICTION_DIAG",
                    "n_left_hits": len(Lhits),
                    "n_right_hits": len(Rhits),
                    "n_pairs_total": n_pairs_total,
                })
            continue

        # Build valid candidates across all chromosomes
        valid_candidates: List[Dict[str, object]] = []
        candidates_by_combo: Dict[Tuple[int, int], List[Dict[str, object]]] = defaultdict(list)

        for lo, ro in expected_combos:
            Lhits = left_hits_by_off.get(lo, [])
            Rhits = right_hits_by_off.get(ro, [])
            n_pairs_total = len(Lhits) * len(Rhits)

            if n_pairs_total == 0:
                # Still useful to emit UNSUPPORTED diagnostics later
                continue

            for L in Lhits:
                for R in Rhits:
                    # Candidate must be on the same chromosome
                    if L.chrom != R.chrom:
                        continue

                    # Strict order consistency (Variant A)
                    if pred_strand == "+":
                        # Must have L before R: L.end <= R.start
                        if L.end > R.start:
                            continue
                        g_start = int(L.end)
                        g_end = int(R.start)
                    else:
                        # For '-' strand, we expect R before L: R.end <= L.start
                        if R.end > L.start:
                            continue
                        g_start = int(R.end)
                        g_end = int(L.start)

                    if g_end <= g_start:
                        continue

                    gap_chr = L.chrom
                    gap_len = max(0, g_end - g_start)

                    # Overlap/recall computed against predicted P only if same chromosome
                    if gap_chr == pred_chr_hal:
                        ov = overlap_len(pred_start, pred_end, g_start, g_end)
                        iou = iou_1d(pred_start, pred_end, g_start, g_end)
                    else:
                        ov = 0
                        iou = 0.0  # By definition, different chromosome => IoU=0

                    pred_recall = (ov / pred_len) if pred_len > 0 else 0.0
                    lr = compute_len_ratio(gap_len, orig_len)

                    row = {
                        "group_id": gid,
                        "pred_chr_hal": pred_chr_hal,
                        "pred_strand": pred_strand,
                        "pred_start": pred_start,
                        "pred_end": pred_end,
                        "pred_len": pred_len,
                        "orig_len": orig_len,
                        "left_offset": lo,
                        "right_offset": ro,
                        "gap_chr": gap_chr,
                        "gap_start": g_start,
                        "gap_end": g_end,
                        "gap_len": gap_len,
                        "overlap_bp": int(ov),
                        "pred_recall": float(pred_recall),
                        "iou": float(iou),
                        "len_ratio": float(lr),
                        "class": "",  # to be filled later
                        "n_left_hits": len(Lhits),
                        "n_right_hits": len(Rhits),
                        "n_pairs_total": n_pairs_total,
                    }
                    valid_candidates.append(row)
                    candidates_by_combo[(lo, ro)].append(row)

        n_valid = len(valid_candidates)

        # UNSUPPORTED: no valid candidates at all (after base validity: same chr + orientation + g_end>g_start)
        if n_valid == 0:
            wm.writerow({
                "group_id": gid,
                "pred_chr_in": pred_chr_in,
                "pred_chr_hal": pred_chr_hal,
                "pred_strand": pred_strand,
                "pred_start": pred_start,
                "pred_end": pred_end,
                "pred_len": pred_len,
                "orig_len": orig_len,
                "hal_status": "HAL_UNSUPPORTED",
                "n_valid_candidates": 0,
                "n_support": 0,
                "n_alt_other_real": 0,
                "best_left_offset": "",
                "best_right_offset": "",
                "best_gap_chr": "",
                "best_gap_start": "",
                "best_gap_end": "",
                "best_gap_len": "",
                "best_overlap_bp": 0,
                "best_pred_recall": 0.0,
                "best_iou": 0.0,
                "best_len_ratio": "",
                "R_min_pred_recall": args.min_pred_recall,
                "I_alt_iou_threshold": args.alt_iou_threshold,
                "len_ratio_max": args.len_ratio_max,
                "note": "NO_VALID_CANDIDATES",
            })

            # Emit per-combo diagnostics into out_unsup (hit counts + pairs_total)
            for lo, ro in expected_combos:
                Lhits = left_hits_by_off.get(lo, [])
                Rhits = right_hits_by_off.get(ro, [])
                n_pairs_total = len(Lhits) * len(Rhits)
                wu.writerow({
                    "group_id": gid,
                    "pred_chr_hal": pred_chr_hal,
                    "pred_strand": pred_strand,
                    "pred_start": pred_start,
                    "pred_end": pred_end,
                    "pred_len": pred_len,
                    "orig_len": orig_len,
                    "left_offset": lo,
                    "right_offset": ro,
                    "gap_chr": "",
                    "gap_start": "",
                    "gap_end": "",
                    "gap_len": "",
                    "overlap_bp": 0,
                    "pred_recall": 0.0,
                    "iou": 0.0,
                    "len_ratio": "",
                    "class": "UNSUPPORTED_DIAG",
                    "n_left_hits": len(Lhits),
                    "n_right_hits": len(Rhits),
                    "n_pairs_total": n_pairs_total,
                })
            continue

        # Classify candidates under the new logic
        R = float(args.min_pred_recall)
        I = float(args.alt_iou_threshold)
        len_ratio_max = float(args.len_ratio_max)

        support: List[Dict[str, object]] = []
        alt_other_real: List[Dict[str, object]] = []
        alt: List[Dict[str, object]] = []

        for c in valid_candidates:
            pr = float(c["pred_recall"])
            if pr >= R:
                c["class"] = "SUPPORT"
                support.append(c)
            else:
                c["class"] = "ALT"
                alt.append(c)

                # AltOther: "geometrically other area" if IoU < I
                if float(c["iou"]) < I:
                    # AltOtherReal: also comparable in length to orig_len
                    if float(c["len_ratio"]) <= len_ratio_max:
                        c["class"] = "ALT_OTHER_REAL"
                        alt_other_real.append(c)
                    else:
                        c["class"] = "ALT_OTHER_LONGMISMATCH"

        n_support = len(support)
        n_alt_other_real = len(alt_other_real)

        # Determine status
        if n_support == 0:
            hal_status = "HAL_DISAGREE"
            note = "SUPPORT_EMPTY"
        else:
            if n_alt_other_real == 0:
                hal_status = "HAL_CONFIRMED_UNIQUE"
                note = "SUPPORT_PRESENT_NO_ALT_OTHER_REAL"
            else:
                hal_status = "HAL_AMBIGUOUS"
                note = f"SUPPORT_PRESENT_AND_ALT_OTHER_REAL(n={n_alt_other_real})"

        # Pick best candidate for reporting:
        # best is chosen among SUPPORT only; if no SUPPORT, best fields remain empty/zero.
        best = None
        if n_support > 0:
            best = choose_best_support(support)

        wm.writerow({
            "group_id": gid,
            "pred_chr_in": pred_chr_in,
            "pred_chr_hal": pred_chr_hal,
            "pred_strand": pred_strand,
            "pred_start": pred_start,
            "pred_end": pred_end,
            "pred_len": pred_len,
            "orig_len": orig_len,
            "hal_status": hal_status,
            "n_valid_candidates": n_valid,
            "n_support": n_support,
            "n_alt_other_real": n_alt_other_real,
            "best_left_offset": int(best["left_offset"]) if best is not None else "",
            "best_right_offset": int(best["right_offset"]) if best is not None else "",
            "best_gap_chr": best["gap_chr"] if best is not None else "",
            "best_gap_start": int(best["gap_start"]) if best is not None else "",
            "best_gap_end": int(best["gap_end"]) if best is not None else "",
            "best_gap_len": int(best["gap_len"]) if best is not None else "",
            "best_overlap_bp": int(best["overlap_bp"]) if best is not None else 0,
            "best_pred_recall": float(best["pred_recall"]) if best is not None else 0.0,
            "best_iou": float(best["iou"]) if best is not None else 0.0,
            "best_len_ratio": float(best["len_ratio"]) if best is not None else "",
            "R_min_pred_recall": R,
            "I_alt_iou_threshold": I,
            "len_ratio_max": len_ratio_max,
            "note": note,
        })

        # Candidate dumps:
        # - AMBIGUOUS: all SUPPORT + all ALT_OTHER_REAL
        # - DISAGREE: all valid candidates with pred_recall < R (ALT and its subclasses)
        # - OTHERCHR: all valid candidates where gap_chr != pred_chr_hal (useful debug; independent of status)
        if hal_status == "HAL_AMBIGUOUS":
            for c in support:
                wa.writerow(c)
            for c in alt_other_real:
                wa.writerow(c)

        if hal_status == "HAL_DISAGREE":
            for c in valid_candidates:
                if float(c["pred_recall"]) < R:
                    wd.writerow(c)

        # Always emit other-chrom candidates into out_other (helps inspect cross-chr ambiguity/disagree)
        for c in valid_candidates:
            if str(c["gap_chr"]) != pred_chr_hal:
                wo.writerow(c)

    out_main.close()
    out_amb.close()
    out_dis.close()
    out_other.close()
    out_unsup.close()


if __name__ == "__main__":
    main()
