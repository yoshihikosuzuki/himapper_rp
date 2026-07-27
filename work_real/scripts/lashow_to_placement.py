#!/usr/bin/env python3
"""Convert himapper LAshow output to a normalized mapping table.

LAshow gives the A-read (mapped read) and B-read (reference contig) by 1-based INDEX
in their respective DAMs.  reads.dam is in chr6.reads.fa order and ref.dam in
chr6panel.fa order, so we map indices back to names via the two ordered name lists.

By default (primary): one placement per read — the longest-span alignment
  columns: read_name  hap  ctg  start  end  strand
With --all: EVERY candidate alignment (primary + secondaries), for secondary analysis
  columns: read_name  hap  ctg  start  end  strand  span

hap = mat/pat from the contig-name suffix (_MATERNAL/_PATERNAL).

Usage:
  LAshow reads.dam ref.dam reads.ref.las | \\
    lashow_to_placement.py [--all] <reads.names> <ctg.names> > mapping.tsv
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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true",
                    help="emit every candidate alignment (+span col), not just the longest")
    ap.add_argument("reads_names", help="read names, one per line, in DAM/fasta order")
    ap.add_argument("ctg_names", help="contig names, one per line, in DAM/fasta order")
    args = ap.parse_args()
    reads = load_names(args.reads_names)
    ctgs = load_names(args.ctg_names)

    best = {}  # read_name -> (span, hap, ctg, start, end, strand)   (primary mode only)
    for line in sys.stdin:
        m = LINE.match(line)
        if not m:
            continue
        a_idx = num(m.group(1)) - 1
        b_idx = num(m.group(3)) - 1
        strand = "+" if m.group(5) == "n" else "-"
        b1, b2 = num(m.group(8)), num(m.group(9))
        s, e = min(b1, b2), max(b1, b2)
        if not (0 <= a_idx < len(reads)) or not (0 <= b_idx < len(ctgs)):
            continue
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
