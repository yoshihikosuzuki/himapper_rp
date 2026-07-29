#!/usr/bin/env python3
"""Build the HG002v1.1 rDNA/GAP + active-HOR evaluation mask.

The rDNA/GAP input is already a set of per-array envelopes.  From the CenSat
input, only records whose feature name starts with ``active_hor`` are kept.
All selected 0-based, half-open intervals are sorted and merged before the BED
and provenance metadata are written atomically.
"""

import argparse
import hashlib
import os
import platform
import shlex
import sys
import tempfile


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_bed(path, *, feature_prefix=None):
    records = []
    with open(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: fewer than 3 columns")
            if feature_prefix is not None:
                if len(fields) < 4 or not fields[3].startswith(feature_prefix):
                    continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: non-integer BED coordinate"
                ) from error
            if start < 0 or end <= start:
                raise ValueError(
                    f"{path}:{line_number}: invalid interval "
                    f"{fields[0]}:{start}-{end}"
                )
            records.append((fields[0], start, end))
    if not records:
        selection = (
            f" with feature prefix {feature_prefix!r}"
            if feature_prefix is not None
            else ""
        )
        raise RuntimeError(f"no BED intervals loaded from {path}{selection}")
    return records


def merge_records(records):
    merged = []
    for contig, start, end in sorted(records):
        if (
            merged
            and merged[-1][0] == contig
            and start <= merged[-1][2]
        ):
            old_contig, old_start, old_end = merged[-1]
            merged[-1] = (old_contig, old_start, max(old_end, end))
        else:
            merged.append((contig, start, end))
    return tuple(merged)


def atomic_write(path, write_callback):
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{os.path.basename(path)}.", dir=directory, text=True
    )
    try:
        with os.fdopen(descriptor, "w") as handle:
            write_callback(handle)
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rdna-gap-bed", required=True)
    parser.add_argument("--censat-bed", required=True)
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--expected-intervals", type=int)
    parser.add_argument("--expected-bp", type=int)
    args = parser.parse_args()

    rdna_records = load_bed(args.rdna_gap_bed)
    active_records = load_bed(args.censat_bed, feature_prefix="active_hor")
    rdna_merged = merge_records(rdna_records)
    active_merged = merge_records(active_records)
    combined = merge_records((*rdna_merged, *active_merged))
    combined_bp = sum(end - start for _contig, start, end in combined)

    if (
        args.expected_intervals is not None
        and len(combined) != args.expected_intervals
    ):
        raise RuntimeError(
            f"combined interval count {len(combined)} != expected "
            f"{args.expected_intervals}"
        )
    if args.expected_bp is not None and combined_bp != args.expected_bp:
        raise RuntimeError(
            f"combined mask bp {combined_bp} != expected {args.expected_bp}"
        )

    atomic_write(
        args.output_bed,
        lambda handle: handle.writelines(
            f"{contig}\t{start}\t{end}\t"
            "rDNA_GAP_or_active_HOR\n"
            for contig, start, end in combined
        ),
    )
    output_sha256 = sha256(args.output_bed)

    metadata_rows = (
        ("command", shlex.join(sys.argv)),
        ("python", platform.python_version()),
        ("coordinate_system", "0-based_half-open"),
        ("reference", "T2T-HG002v1.1"),
        ("rdna_gap_source", os.path.realpath(args.rdna_gap_bed)),
        ("rdna_gap_source_sha256", sha256(args.rdna_gap_bed)),
        ("rdna_gap_source_records", str(len(rdna_records))),
        ("rdna_gap_merged_intervals", str(len(rdna_merged))),
        (
            "rdna_gap_bp",
            str(sum(end - start for _contig, start, end in rdna_merged)),
        ),
        ("censat_source", os.path.realpath(args.censat_bed)),
        ("censat_source_sha256", sha256(args.censat_bed)),
        ("censat_filter", "column4 startswith active_hor"),
        ("active_hor_source_records", str(len(active_records))),
        ("active_hor_merged_intervals", str(len(active_merged))),
        (
            "active_hor_bp",
            str(sum(end - start for _contig, start, end in active_merged)),
        ),
        ("merge_rule", "same-contig start<=current_end"),
        ("combined_intervals", str(len(combined))),
        ("combined_bp", str(combined_bp)),
        ("output_bed", os.path.realpath(args.output_bed)),
        ("output_bed_sha256", output_sha256),
        (
            "evaluation_exclusion_rule",
            "any_1bp_truth_interval_overlap",
        ),
    )
    atomic_write(
        args.metadata,
        lambda handle: handle.writelines(
            f"{key}\t{value}\n" for key, value in metadata_rows
        ),
    )

    print(
        f"wrote {args.output_bed}: {len(combined)} intervals, "
        f"{combined_bp} bp"
    )


if __name__ == "__main__":
    main()
