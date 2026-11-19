#!/usr/bin/env python3
"""Export selected RepeatMasker annotations, or all external intervals, as BED6."""

import argparse
import json
from annotation_io import BedRepeat, iter_bed, write_bed6_row


def repeatmasker_out_to_bed(input_path, output_path, target_classes=("LTR",),
                          target_repeats=(), chrom_map=None):
    """Preserve the notebook converter's class-prefix and coordinate rules."""
    if isinstance(target_classes, str):
        target_classes = [target_classes]
    target_classes = tuple(target_classes)
    with open(input_path) as infile, open(output_path, "w") as outfile:
        for line in infile:
            fields = line.split()
            if not fields or not fields[0].isdigit() or len(fields) < 15:
                continue
            rep_class, name = fields[10], fields[9]
            keep = ("all" in target_classes or name in target_repeats or
                    any(rep_class == target or rep_class.startswith(target + "/")
                        for target in target_classes))
            if not keep:
                continue
            if fields[8] not in {"+", "C"}:
                raise ValueError(f"Unexpected strand value: {fields[8]}")
            row = BedRepeat(fields[4], int(fields[5]) - 1, int(fields[6]),
                            "+" if fields[8] == "+" else "-", name)
            write_bed6_row(outfile, row, chrom_map)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--format", choices=["rm_out", "bed6", "custom5"], required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--target-groups", default="LTR")
    parser.add_argument("--target-repeats", default="")
    parser.add_argument("--chrom-map", default="")
    args = parser.parse_args()
    chrom_map = {}
    if args.chrom_map:
        with open(args.chrom_map) as handle:
            chrom_map = json.load(handle)
    if args.format == "rm_out":
        groups = [x.strip() for x in args.target_groups.split(",") if x.strip()]
        names = {x.strip() for x in args.target_repeats.split(",") if x.strip()}
        if not groups and not names:
            parser.error("At least one target group or repeat name is required for rm_out")
        repeatmasker_out_to_bed(args.input, args.output, groups, names, chrom_map)
    else:
        with open(args.output, "w") as handle:
            for row in iter_bed(args.input, args.format):
                write_bed6_row(handle, row, chrom_map)


if __name__ == "__main__":
    main()
