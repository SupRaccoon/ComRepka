#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Prepare ComRepka-like flanking queries around target repeats.

What this script does:
  1. Reads target repeat coordinates from RepeatMasker .out.
  2. For every target repeat, builds one genomic-left/upstream query and one genomic-right/downstream query.
  3. The left/upstream query follows the ComRepka non-repeat tag idea:
       - look inside the upstream window for the nearest non-repeat interval >= --nrt-len,
       - if found, start the upstream query from that interval,
       - include everything from that interval to the target repeat,
       - add --overhang-len bp inside the target repeat.
     If no such non-repeat interval is found, use the full upstream window.
  4. The right/downstream query is the full downstream window plus --overhang-len bp inside the target repeat.
  5. Writes output files in the same column/header style as the provided flank-preparation script.

Coordinate convention in TSV/FASTA headers:
  - 0-based, half-open intervals: [start, end).
  - left_flank_start/end and right_flank_start/end describe the actual FASTA query sequence,
    including the repeat-overhang segment.

Notes:
  - "upstream" and "downstream" here mean genomic-left and genomic-right, matching the previous script's
    left/right output format. They are not reversed for repeats annotated on the minus strand.
  - Soft masking is optional. If enabled, RepeatMasker intervals inside the query are lower-cased, except
    repeat ends, target-neighbourhood windows, and overhangs, which remain upper-case.
"""

import argparse
import json
import os
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
import pandas as pd
from annotation_io import iter_bed
from tqdm import tqdm


# -----------------------------
# FASTA random access via .fai
# -----------------------------
@dataclass
class FaiEntry:
    length: int
    offset: int
    line_bases: int
    line_width: int


def load_fai(fai_path: str) -> Dict[str, FaiEntry]:
    fai: Dict[str, FaiEntry] = {}
    with open(fai_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            name, ln, off, lb, lw = line.rstrip("\n").split("\t")[:5]
            fai[name] = FaiEntry(int(ln), int(off), int(lb), int(lw))
    return fai


def faidx_fetch(fasta_path: str, fai: Dict[str, FaiEntry], chrom: str, start0: int, end0: int) -> str:
    """
    Fetch sequence from FASTA using .fai index.
    Coordinates are 0-based half-open [start0, end0).
    Returns uppercase sequence without newlines.
    """
    if chrom not in fai:
        raise KeyError(f"Chromosome '{chrom}' not found in FASTA index")
    e = fai[chrom]
    if start0 < 0 or end0 < 0 or end0 < start0:
        raise ValueError(f"Bad coords: {chrom}:{start0}-{end0}")
    if end0 > e.length:
        raise ValueError(f"Out of bounds: {chrom}:{start0}-{end0}, chr_len={e.length}")

    def bytepos(pos0: int) -> int:
        line_no = pos0 // e.line_bases
        in_line = pos0 % e.line_bases
        return e.offset + line_no * e.line_width + in_line

    start_byte = bytepos(start0)
    end_byte = bytepos(end0)
    need = end0 - start0

    with open(fasta_path, "rb") as fh:
        fh.seek(start_byte)
        raw = fh.read(end_byte - start_byte)

    seq = raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()
    if len(seq) < need:
        with open(fasta_path, "rb") as fh:
            fh.seek(start_byte)
            raw = fh.read((end_byte - start_byte) + (need - len(seq)) + 1024)
        seq = raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()
    return seq[:need]


# -----------------------------
# RepeatMasker .out parsing
# -----------------------------
@dataclass
class RMRow:
    chrom: str
    start0: int
    end0: int
    strand: str
    rep_name: str
    class_family: str
    rm_id: str


def iter_rm_out(path: str) -> Iterable[RMRow]:
    """
    Parse RepeatMasker .out file.
    RepeatMasker begin/end are 1-based inclusive.
    Converted to 0-based half-open: start0=begin-1, end0=end.
    """
    with open(path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            s = line.lstrip()
            if s.startswith("bit") or s.startswith("score"):
                continue
            if not s or not s[0].isdigit():
                continue

            parts = line.split()
            if len(parts) < 15:
                continue

            rm_id: Optional[str] = None
            for tok in reversed(parts):
                if tok == "*":
                    continue
                if tok.isdigit():
                    rm_id = tok
                    break
            if rm_id is None:
                continue

            chrom = parts[4]
            begin = int(parts[5])
            end = int(parts[6])
            strand_raw = parts[8]
            strand = "-" if strand_raw in {"C", "-"} else "+"
            rep_name = parts[9]
            class_family = parts[10]

            start0 = begin - 1
            end0 = end
            if end0 <= start0:
                continue

            yield RMRow(
                chrom=chrom,
                start0=start0,
                end0=end0,
                strand=strand,
                rep_name=rep_name,
                class_family=class_family,
                rm_id=rm_id,
            )


# -----------------------------
# Interval helpers
# -----------------------------
def merge_intervals(iv: Sequence[Tuple[int, int]], join_gap: int = 0) -> List[Tuple[int, int]]:
    """Merge sorted or unsorted 0-based half-open intervals."""
    if not iv:
        return []
    items = sorted(iv, key=lambda x: (x[0], x[1]))
    out: List[Tuple[int, int]] = []
    cs, ce = items[0]
    for s, e in items[1:]:
        if s <= ce + join_gap:
            ce = max(ce, e)
        else:
            out.append((cs, ce))
            cs, ce = s, e
    out.append((cs, ce))
    return out


def clipped_overlaps(intervals: Sequence[Tuple[int, int]], start0: int, end0: int) -> List[Tuple[int, int]]:
    out = []
    for s, e in intervals:
        if e <= start0:
            continue
        if s >= end0:
            break
        cs = max(s, start0)
        ce = min(e, end0)
        if ce > cs:
            out.append((cs, ce))
    return out


def complement_gaps(start0: int, end0: int, blocked: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Return non-blocked intervals inside [start0, end0)."""
    gaps: List[Tuple[int, int]] = []
    cur = start0
    for s, e in blocked:
        if s > cur:
            gaps.append((cur, s))
        cur = max(cur, e)
    if cur < end0:
        gaps.append((cur, end0))
    return gaps


def find_upstream_query_start(
    repeat_start0: int,
    upstream_window_start0: int,
    repeat_intervals_chr: Sequence[Tuple[int, int]],
    nrt_len: int,
) -> Tuple[int, str, int]:
    """
    Find ComRepka-like NRT start for upstream query.

    Returns:
      query_start0, nrt_status, nrt_len_observed

    nrt_status:
      - immediate_nrt: nearest non-repeat interval >= nrt_len touches the target repeat start
      - distant_nrt: nearest non-repeat interval >= nrt_len exists but is separated from target by repeats/small gaps
      - no_nrt: no suitable non-repeat interval found, full upstream window is used
    """
    if repeat_start0 <= upstream_window_start0:
        return upstream_window_start0, "no_upstream_window", 0

    overlaps = clipped_overlaps(repeat_intervals_chr, upstream_window_start0, repeat_start0)
    blocked = merge_intervals(overlaps, join_gap=0)
    gaps = complement_gaps(upstream_window_start0, repeat_start0, blocked)

    # Scan from the target-proximal side towards the 5' side of the upstream flank.
    for gs, ge in reversed(gaps):
        gap_len = ge - gs
        if gap_len >= nrt_len:
            status = "immediate_nrt" if ge == repeat_start0 else "distant_nrt"
            return gs, status, gap_len

    return upstream_window_start0, "no_nrt", 0


def apply_softmask(
    seq: str,
    query_start0: int,
    query_end0: int,
    repeat_intervals_chr: Sequence[Tuple[int, int]],
    target_start0: int,
    target_end0: int,
    repeat_end_keep: int,
    target_near_keep: int,
) -> str:
    """
    Lowercase RepeatMasker intervals inside the query, then uppercase protected regions.
    Repeat ends, target neighbourhoods and target overhangs remain uppercase.
    """
    chars = list(seq.upper())

    def set_case(abs_s: int, abs_e: int, upper: bool) -> None:
        s = max(abs_s, query_start0)
        e = min(abs_e, query_end0)
        if e <= s:
            return
        rs = s - query_start0
        re = e - query_start0
        if upper:
            chars[rs:re] = [c.upper() for c in chars[rs:re]]
        else:
            chars[rs:re] = [c.lower() for c in chars[rs:re]]

    # 1) Lowercase all repeat intervals.
    for s, e in clipped_overlaps(repeat_intervals_chr, query_start0, query_end0):
        set_case(s, e, upper=False)

    # 2) Uppercase 5' and 3' ends of every repeat.
    for s, e in clipped_overlaps(repeat_intervals_chr, query_start0, query_end0):
        set_case(s, min(e, s + repeat_end_keep), upper=True)
        set_case(max(s, e - repeat_end_keep), e, upper=True)

    # 3) Uppercase target-neighbourhood windows.
    set_case(target_start0 - target_near_keep, target_start0, upper=True)
    set_case(target_end0, target_end0 + target_near_keep, upper=True)

    # 4) Uppercase the target itself within the query overhangs.
    set_case(target_start0, target_end0, upper=True)

    return "".join(chars)


# -----------------------------
# Target selection
# -----------------------------
def split_csv(value: str) -> set:
    return {x.strip() for x in value.split(",") if x.strip()}


def is_target(row: RMRow, target_groups: set, target_repeats: set) -> bool:
    """
    Target definition:
      - repeat name exact match in --target-repeats, OR
      - class/family exact match in --target-groups, OR
      - class prefix before '/' match in --target-groups.
    """
    if target_repeats and row.rep_name in target_repeats:
        return True
    if target_groups:
        if row.class_family in target_groups:
            return True
        prefix = row.class_family.split("/")[0]
        if prefix in target_groups:
            return True
    return False


def wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Prepare ComRepka-like NRT upstream and full downstream flanking queries around target repeats."
    )
    inputs = ap.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--rm-out", help="RepeatMasker .out file for source/main genome")
    inputs.add_argument("--bed", help="External repeat intervals; all rows are targets and masking intervals")
    ap.add_argument("--bed-format", choices=["bed6", "custom5"], default="bed6")
    ap.add_argument("--genome-fasta", required=True, help="Source/main genome FASTA, indexed with samtools faidx")

    ap.add_argument("--target-groups", default="LTR",
                    help="Comma-separated target class/family list, e.g. 'LTR,LTR/ERV1'. Optional if --target-repeats is used.")
    ap.add_argument("--target-repeats", default="",
                    help="Comma-separated exact RepeatMasker repeat names, e.g. 'L1HS,AluYa5'. Optional if --target-groups is used.")

    ap.add_argument("--flank-len", type=int, default=5000,
                    help="Maximum flanking window length outside the target repeat, bp; ComRepka default is 5000")
    ap.add_argument("--nrt-len", type=int, default=500,
                    help="Minimum non-repeat tag length for upstream query, bp; ComRepka default is 500")
    ap.add_argument("--overhang-len", type=int, default=50,
                    help="Repeat-overhang length added inside target repeat on both sides, bp; ComRepka default is 50")
    ap.add_argument("--soft-mask", action="store_true",
                    help="Soft-mask repeat intervals in output FASTA, approximating ComRepka masking")
    ap.add_argument("--target-near-keep", type=int, default=500,
                    help="When --soft-mask is used, keep this many bp near target in uppercase")
    ap.add_argument("--repeat-end-keep", type=int, default=50,
                    help="When --soft-mask is used, keep this many bp at each repeat end in uppercase")

    ap.add_argument("--max-groups", type=int, default=0, help="Debug: limit target repeats written; 0 = no limit")
    ap.add_argument("--out-table", default="flanks.tsv", help="Output TSV with flank table")
    ap.add_argument("--out-fasta", default="flanks.fa", help="Output FASTA with flank query sequences")
    ap.add_argument("--out-targets", default="targets_in_groups.tsv", help="Output TSV: target repeats")
    ap.add_argument("--chrom-map", default="", help="Optional JSON mapping for rm-out chrom -> FASTA chrom")

    args = ap.parse_args()

    target_groups = split_csv(args.target_groups)
    target_repeats = split_csv(args.target_repeats)
    if args.rm_out and not target_groups and not target_repeats:
        raise SystemExit("At least one of --target-groups or --target-repeats must be non-empty")

    if args.rm_out and args.target_groups == 'all':
        target_groups = set(pd.read_csv(args.rm_out, skiprows=2, sep=r'\s+', usecols=[10], names=['rep_name'], dtype=str).rep_name.unique())
    
    if args.flank_len < 0 or args.nrt_len < 1 or args.overhang_len < 0:
        raise SystemExit("Bad lengths: require --flank-len >= 0, --nrt-len >= 1, --overhang-len >= 0")
    print(target_groups)
    chrom_map = {}
    if args.chrom_map:
        with open(args.chrom_map, "r") as f:
            chrom_map = json.load(f)

    fasta = args.genome_fasta
    fai_path = fasta + ".fai"
    if not os.path.exists(fai_path):
        raise SystemExit(f"FASTA index not found: {fai_path}. Create it with: samtools faidx {fasta}")
    fai = load_fai(fai_path)

    repeats_by_chr: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    targets_by_chr: Dict[str, List[RMRow]] = defaultdict(list)

    print("[START] Parse repeat annotation")
    n_rows = 0
    n_kept_rows = 0
    n_targets = 0
    rows = iter_rm_out(args.rm_out) if args.rm_out else iter_bed(args.bed, args.bed_format)
    seen_targets = set()
    for row in tqdm(rows, unit="rows"):
        n_rows += 1
        chrom = chrom_map.get(row.chrom, row.chrom)
        if chrom not in fai:
            continue
        mapped = RMRow(
            chrom=chrom,
            start0=row.start0,
            end0=row.end0,
            strand=row.strand,
            rep_name=row.rep_name,
            class_family=row.class_family,
            rm_id=row.rm_id,
        )
        repeats_by_chr[chrom].append((mapped.start0, mapped.end0))
        n_kept_rows += 1
        if args.bed or is_target(mapped, target_groups, target_repeats):
            key = (chrom, mapped.start0, mapped.end0)
            if key in seen_targets:
                raise ValueError(f"Duplicate target coordinates would produce duplicate group_id: {key}")
            seen_targets.add(key)
            if not 0 <= mapped.start0 < mapped.end0 <= fai[chrom].length:
                raise ValueError(f"Target outside FASTA bounds: {key}")
            targets_by_chr[chrom].append(mapped)
            n_targets += 1

    for chrom in list(repeats_by_chr.keys()):
        repeats_by_chr[chrom] = merge_intervals(repeats_by_chr[chrom], join_gap=0)
        targets_by_chr[chrom].sort(key=lambda r: (r.start0, r.end0, r.rep_name, r.rm_id))

    print(f"[DONE ] Parsed rows={n_rows}, rows_in_fasta={n_kept_rows}, targets={n_targets}, chromosomes={len(repeats_by_chr)}")

    table_header = [
        "group_id",
        "chrom",
        "group_start",
        "group_end",
        "blocked_start",
        "blocked_end",
        "left_flank_start",
        "left_flank_end",
        "right_flank_start",
        "right_flank_end",
        "n_repeats_inside",
        "repeat_ids_inside",
    ]
    targets_header = [
        "group_id", "chrom", "repeat_start", "repeat_end", "repeat_strand",
        "repeat_name", "class_family", "rm_id"
    ]

    print("[START] Build ComRepka-like flanking queries")
    group_count = 0

    with open(args.out_table, "w") as tsv, open(args.out_targets, "w") as tsv2, open(args.out_fasta, "w") as fa:
        tsv.write("\t".join(table_header) + "\n")
        tsv2.write("\t".join(targets_header) + "\n")

        for chrom in tqdm(sorted(targets_by_chr.keys()), unit="chr"):
            chr_len = fai[chrom].length
            repeat_intervals = repeats_by_chr[chrom]

            for target in tqdm(targets_by_chr[chrom]):
                repeat_len = target.end0 - target.start0
                overhang = min(args.overhang_len, repeat_len)

                upstream_window_start = max(0, target.start0 - args.flank_len)
                upstream_query_start, nrt_status, nrt_observed_len = find_upstream_query_start(
                    repeat_start0=target.start0,
                    upstream_window_start0=upstream_window_start,
                    repeat_intervals_chr=repeat_intervals,
                    nrt_len=args.nrt_len,
                )
                upstream_query_end = min(chr_len, target.start0 + overhang)

                downstream_query_start = max(0, target.end0 - overhang)
                downstream_query_end = min(chr_len, target.end0 + args.flank_len)

                group_id = f"{chrom}:{target.start0}-{target.end0}"

                # targets TSV: same column format as the sample; one target per group.
                tsv2.write("\t".join([
                    group_id,
                    chrom,
                    str(target.start0),
                    str(target.end0),
                    target.strand,
                    target.rep_name,
                    target.class_family,
                    target.rm_id,
                ]) + "\n")

                row = {
                    "group_id": group_id,
                    "chrom": chrom,
                    "group_start": str(target.start0),
                    "group_end": str(target.end0),
                    "blocked_start": str(target.start0),
                    "blocked_end": str(target.end0),
                    "left_flank_start": str(upstream_query_start),
                    "left_flank_end": str(upstream_query_end),
                    "right_flank_start": str(downstream_query_start),
                    "right_flank_end": str(downstream_query_end),
                    "n_repeats_inside": "1",
                    "repeat_ids_inside": target.rep_name,
                }
                tsv.write("\t".join(row[h] for h in table_header) + "\n")

                seq_left = faidx_fetch(fasta, fai, chrom, upstream_query_start, upstream_query_end)
                seq_right = faidx_fetch(fasta, fai, chrom, downstream_query_start, downstream_query_end)

                if args.soft_mask:
                    seq_left = apply_softmask(
                        seq=seq_left,
                        query_start0=upstream_query_start,
                        query_end0=upstream_query_end,
                        repeat_intervals_chr=repeat_intervals,
                        target_start0=target.start0,
                        target_end0=target.end0,
                        repeat_end_keep=args.repeat_end_keep,
                        target_near_keep=args.target_near_keep,
                    )
                    seq_right = apply_softmask(
                        seq=seq_right,
                        query_start0=downstream_query_start,
                        query_end0=downstream_query_end,
                        repeat_intervals_chr=repeat_intervals,
                        target_start0=target.start0,
                        target_end0=target.end0,
                        repeat_end_keep=args.repeat_end_keep,
                        target_near_keep=args.target_near_keep,
                    )

                # FASTA header keeps the same field style as the sample script.
                # Extra fields are appended after block:..., so parsers that split the first fields still work.
                left_header = (
                    f">{group_id}|{chrom}|left|L={len(seq_left)}|left:{upstream_query_start}-{upstream_query_end}"
                    f"|block:{target.start0}-{target.end0}|nrt:{nrt_status}|nrt_len:{nrt_observed_len}"
                    f"|overhang:{overhang}"
                )
                right_header = (
                    f">{group_id}|{chrom}|right|L={len(seq_right)}|right:{downstream_query_start}-{downstream_query_end}"
                    f"|block:{target.start0}-{target.end0}|overhang:{overhang}"
                )

                fa.write(left_header + "\n" + wrap_fasta(seq_left) + "\n")
                fa.write(right_header + "\n" + wrap_fasta(seq_right) + "\n")

                group_count += 1
                if args.max_groups and group_count >= args.max_groups:
                    break

            if args.max_groups and group_count >= args.max_groups:
                break

    print(f"[DONE ] Wrote groups={group_count}")
    print(f"[OK] table:  {args.out_table}")
    print(f"[OK] fasta:  {args.out_fasta}")
    print(f"[OK] targets:{args.out_targets}")


if __name__ == "__main__":
    main()
