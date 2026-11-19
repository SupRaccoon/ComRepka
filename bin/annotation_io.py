"""Read explicitly selected repeat interval formats without filtering records."""

from dataclasses import dataclass


@dataclass(frozen=True)
class BedRepeat:
    chrom: str
    start0: int
    end0: int
    strand: str
    rep_name: str
    class_family: str = ""
    rm_id: str = ""


def iter_bed(path, format_name="bed6"):
    """Yield BED6 or legacy chrom/start/end/strand/name rows; retain duplicates."""
    if format_name not in {"bed6", "custom5"}:
        raise ValueError(f"Unsupported annotation format: {format_name}")
    width = 6 if format_name == "bed6" else 5
    with open(path, encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != width:
                raise ValueError(f"{path}:{line_no}: expected exactly {width} tab-separated fields")
            chrom, start, end = fields[:3]
            if format_name == "bed6":
                name, score, strand = fields[3:]
                if not score.isdigit() or not 0 <= int(score) <= 1000:
                    raise ValueError(f"{path}:{line_no}: BED score must be an integer in 0..1000")
            else:
                strand, name = fields[3:]
            start, end = int(start), int(end)
            if not chrom or any(x.isspace() for x in chrom) or "|" in chrom:
                raise ValueError(f"{path}:{line_no}: invalid chromosome identifier")
            if not name or strand not in {"+", "-", "."} or not 0 <= start < end:
                raise ValueError(f"{path}:{line_no}: invalid name, strand or interval")
            yield BedRepeat(chrom, start, end, strand, name)


def write_bed6_row(handle, row, chrom_map=None):
    chrom = (chrom_map or {}).get(row.chrom, row.chrom)
    handle.write(f"{chrom}\t{row.start0}\t{row.end0}\t{row.rep_name}\t0\t{row.strand}\n")
