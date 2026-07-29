#!/usr/bin/env python3
"""Convert himapper LAshow output to original FASTA scaffold coordinates.

LAshow identifies DAM sequences as 1-based ``scaffold.contig`` pairs.  A FASTA
scaffold containing N runs is represented by multiple sequences in a DAM.
Forward LAshow coordinates already include the sequence's ``fpulse`` offset.
For reverse hits, however, LAshow assumes that the whole scaffold was reverse
complemented, whereas himapper reverse complements each N-delimited DAM contig.
``--ref-layout`` supplies the metadata needed to correct that difference.

By default (primary): one placement per read — the longest-span alignment
  columns: read_name  hap  ctg  start  end  strand
With --all: EVERY candidate alignment (primary + secondaries), for secondary analysis
  columns: read_name  hap  ctg  start  end  strand  span

hap = mat/pat from the contig-name suffix (_MATERNAL/_PATERNAL).

Usage:
  DBshow -n ref.dam > ref.dam.layout
  LAshow reads.dam ref.dam reads.ref.las | \\
    lashow_to_placement.py --ref-layout ref.dam.layout [--all] \\
      <reads.names> <ctg.names> > mapping.tsv
  (reads.names / ctg.names: one name per line, in DAM/fasta order)
"""
import argparse
import re
import sys

NUM = r"[\d,]+"
LINE = re.compile(
    r"^[>+\- ]+"
    r"(" + NUM + r")\.(\d+)\s+"           # A-id . A-ctg   (read index)
    r"(" + NUM + r")\.(\d+)\s+"           # B-id . B-ctg   (contig index)
    r"([nc])\s+"                            # n=fwd, c=rev
    r"[<\[]\s*(" + NUM + r")\s*\.\.\s*"
    r"(" + NUM + r")\s*[>\]]\s+x\s+"        # A coords (ignored)
    r"[<\[]\s*(" + NUM + r")\s*\.\.\s*"
    r"(" + NUM + r")\s*[>\]]"               # B coords (ref start..end)
)
LAYOUT = re.compile(
    r"^>(.*?)\s+::\s+Contig\s+(\d+)"
    r"\[\s*(" + NUM + r")\s*,\s*(" + NUM + r")\s*\]\s*$"
)


def num(s):
    return int(s.replace(",", ""))


def hap_of(ctg):
    if ctg.endswith("_MATERNAL"):
        return "mat"
    if ctg.endswith("_PATERNAL"):
        return "pat"
    return "NA"


def load_names(path):
    with open(path) as fh:
        return [ln.strip().split()[0] for ln in fh if ln.strip()]


def load_ref_layout(path, scaffolds):
    """Load ``DBshow -n`` output keyed by zero-based scaffold/contig IDs."""
    scaffold_ids = {}
    for index, name in enumerate(scaffolds):
        if name in scaffold_ids:
            raise ValueError(f"duplicate scaffold name in name list: {name}")
        scaffold_ids[name] = index

    contigs = {}
    scaffold_lengths = {}
    with open(path) as fh:
        for line_number, line in enumerate(fh, 1):
            match = LAYOUT.match(line.rstrip("\n"))
            if match is None:
                if line.strip():
                    raise ValueError(
                        f"{path}:{line_number}: invalid DBshow layout line"
                    )
                continue
            scaffold_name = match.group(1).split()[0]
            if scaffold_name not in scaffold_ids:
                raise ValueError(
                    f"{path}:{line_number}: unknown scaffold {scaffold_name}"
                )
            scaffold_index = scaffold_ids[scaffold_name]
            contig_index = int(match.group(2))
            fpulse = num(match.group(3))
            end = num(match.group(4))
            if fpulse < 0 or end < fpulse:
                raise ValueError(
                    f"{path}:{line_number}: invalid contig interval "
                    f"[{fpulse},{end}]"
                )
            key = (scaffold_index, contig_index)
            if key in contigs:
                raise ValueError(
                    f"{path}:{line_number}: duplicate scaffold/contig {key}"
                )
            contigs[key] = (fpulse, end - fpulse)
            scaffold_lengths[scaffold_index] = max(
                scaffold_lengths.get(scaffold_index, 0), end
            )

    missing = [
        scaffolds[index]
        for index in range(len(scaffolds))
        if index not in scaffold_lengths
    ]
    if missing:
        raise ValueError(
            f"{path}: no DAM layout records for scaffold(s): "
            + ", ".join(missing)
        )
    for scaffold_index, scaffold in enumerate(scaffolds):
        records = sorted(
            (contig_index, fpulse, contig_length)
            for (index, contig_index), (fpulse, contig_length) in contigs.items()
            if index == scaffold_index
        )
        indices = [record[0] for record in records]
        if indices != list(range(len(records))):
            raise ValueError(
                f"{path}: non-contiguous DAM contig indices for {scaffold}: "
                f"{indices}"
            )
        if records[0][1] != 0:
            raise ValueError(
                f"{path}: unsupported leading-N/all-N scaffold {scaffold}: "
                f"contig 0 starts at {records[0][1]}"
            )
        previous_end = 0
        for contig_index, fpulse, contig_length in records:
            if fpulse < previous_end:
                raise ValueError(
                    f"{path}: overlapping/out-of-order DAM contig "
                    f"{scaffold}.{contig_index}"
                )
            if contig_length == 0 and contig_index != records[-1][0]:
                raise ValueError(
                    f"{path}: zero-length non-final DAM contig "
                    f"{scaffold}.{contig_index}"
                )
            previous_end = fpulse + contig_length
    return contigs, scaffold_lengths


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true",
                    help="emit every candidate alignment (+span col), not just the longest")
    ap.add_argument(
        "--ref-layout",
        required=True,
        help="output of 'DBshow -n ref.dam'",
    )
    ap.add_argument("reads_names", help="read names, one per line, in DAM/fasta order")
    ap.add_argument("ctg_names", help="contig names, one per line, in DAM/fasta order")
    args = ap.parse_args()
    reads = load_names(args.reads_names)
    ctgs = load_names(args.ctg_names)
    contigs, scaffold_lengths = load_ref_layout(args.ref_layout, ctgs)

    best = {}  # read_name -> (span, hap, ctg, start, end, strand)   (primary mode only)
    for line_number, line in enumerate(sys.stdin, 1):
        m = LINE.match(line)
        if not m:
            continue
        a_idx = num(m.group(1)) - 1
        b_idx = num(m.group(3)) - 1
        b_contig = int(m.group(4)) - 1
        strand = "+" if m.group(5) == "n" else "-"
        b1, b2 = num(m.group(8)), num(m.group(9))
        s, e = min(b1, b2), max(b1, b2)
        if not (0 <= a_idx < len(reads)) or not (0 <= b_idx < len(ctgs)):
            continue
        key = (b_idx, b_contig)
        if key not in contigs:
            raise ValueError(
                f"LAshow line {line_number}: missing DAM layout for "
                f"{b_idx + 1}.{b_contig + 1}"
            )
        fpulse, contig_length = contigs[key]
        if strand == "-":
            delta = 2 * fpulse + contig_length - scaffold_lengths[b_idx]
            s += delta
            e += delta
        contig_end = fpulse + contig_length
        if s < fpulse or e > contig_end or e < s:
            raise ValueError(
                f"LAshow line {line_number}: corrected interval [{s},{e}] "
                f"is outside DAM contig {b_idx + 1}.{b_contig + 1} "
                f"interval [{fpulse},{contig_end}]"
            )
        rname, ctg = reads[a_idx], ctgs[b_idx]
        if args.all:
            print(f"{rname}\t{hap_of(ctg)}\t{ctg}\t{s}\t{e}\t{strand}\t{e - s}")
        else:
            span = e - s
            if rname not in best or span > best[rname][0]:
                best[rname] = (span, hap_of(ctg), ctg, s, e, strand)
    if not args.all:
        for rname, (_span, hap, ctg, s, e, strand) in best.items():
            print(f"{rname}\t{hap}\t{ctg}\t{s}\t{e}\t{strand}")


if __name__ == "__main__":
    main()
