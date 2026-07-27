#!/usr/bin/env python3
"""Classify reads by variants in their mapping-truth intervals.

Output: <read_name>\t<no_variant|variant>\t<n_variants>

Usage: variant_strata.py variants.tsv mapping_truth.tsv out.tsv
"""
import bisect
import sys


def main():
    variants = {}
    with open(sys.argv[1]) as fh:
        for line in fh:
            contig, pos = line.rstrip("\n").split("\t")[:2]
            variants.setdefault(contig, []).append(int(pos))
    for positions in variants.values():
        positions.sort()

    with open(sys.argv[2]) as mapping_truth, open(sys.argv[3], "w") as out:
        for line in mapping_truth:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            name, contig = fields[0], fields[2]
            start, end = int(fields[3]), int(fields[4])
            positions = variants.get(contig, ())
            i = bisect.bisect_left(positions, start)
            j = bisect.bisect_left(positions, end)
            n_variants = j - i
            stratum = "variant" if n_variants else "no_variant"
            out.write(f"{name}\t{stratum}\t{n_variants}\n")


if __name__ == "__main__":
    main()
