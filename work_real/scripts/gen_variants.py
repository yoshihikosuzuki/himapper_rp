#!/usr/bin/env python3
"""Extract MAT/PAT haplotype-distinguishing variant positions from a cs-tagged PAF.

Reads one or more `minimap2 -c --cs=long` PAFs (haplotype-vs-haplotype alignment)
and emits every variant (substitution / insertion / deletion) as a position on the
TARGET contig of that PAF.  Run it once per direction so both haplotype contigs get
their own coordinates:
    minimap2 -c --cs=long -x asm5 matctg.fa patctg.fa  -> variants on chr6_MATERNAL
    minimap2 -c --cs=long -x asm5 patctg.fa matctg.fa  -> variants on chr6_PATERNAL

Output (TSV, sorted): <contig>\\t<pos>   (0-based target coordinate of the variant)

Usage: gen_variants.py out.tsv in1.paf [in2.paf ...]
"""
import re
import sys

TOK = re.compile(r'([=:\*\+\-])([A-Za-z0-9]+)')


def variants_from_paf(paf, out_by_ctg):
    for ln in open(paf):
        x = ln.rstrip("\n").split("\t")
        if len(x) < 12:
            continue
        tname, tstart = x[5], int(x[7])
        cs = ""
        for t in x[12:]:
            if t.startswith("cs:Z:"):
                cs = t[5:]
        if not cs:
            continue
        tp = tstart
        acc = out_by_ctg.setdefault(tname, [])
        for m in TOK.finditer(cs):
            op, val = m.group(1), m.group(2)
            if op == "=":
                tp += len(val)
            elif op == ":":
                tp += int(val)
            elif op == "*":            # substitution: 1 target base
                acc.append(tp); tp += 1
            elif op == "+":            # insertion in query (no target advance)
                acc.append(tp)
            elif op == "-":            # deletion from query: spans target bases
                acc.append(tp); tp += len(val)


def main():
    out = sys.argv[1]
    by_ctg = {}
    for paf in sys.argv[2:]:
        variants_from_paf(paf, by_ctg)
    n = 0
    with open(out, "w") as fh:
        for ctg in sorted(by_ctg):
            for pos in sorted(by_ctg[ctg]):
                fh.write(f"{ctg}\t{pos}\n"); n += 1
    print(f"gen_variants: wrote {n} variants over {len(by_ctg)} contigs -> {out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
