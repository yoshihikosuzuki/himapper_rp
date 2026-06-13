#!/usr/bin/env python3
"""Precompute STRICT truth-based resolvability for simulated reads.

A read is *resolvable* iff its true origin can, in principle, be told apart from
every homolog/repeat in the reference -- i.e. its error-free ground-truth
reference span is NOT identical to any other reference locus.  We extract each
read's GT reference span (no sequencing errors), self-map it against the
reference with minimap2, and count loci that align at NM=0 over >=99% of the
span.  Exactly one such locus (the GT itself) => resolvable; >=2 (an identical
homolog exists) => unresolvable.

This is the rigorous counterpart of the lightweight "copy-1 k-mer in the span"
proxy in eval_mapping.py: it also credits spans that are distinguished only by a
copy>=2 (but rearranged) k-mer, which the proxy misses.

Output: <ref_prefix>.strictres.k<K>.txt  with lines "<read_id> <0|1>".

Usage: strict_resolvable.py <ref.fasta> <reads.fasta> [K]
"""
import collections
import os
import re
import subprocess
import sys

K = 32


def parse_ref_seqs(ref_fasta):
    seqs = {}
    idx = 0
    cur = []
    with open(ref_fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs[idx] = "".join(cur)
                    cur = []
                idx += 1
            else:
                cur.append(line.strip())
        if cur:
            seqs[idx] = "".join(cur)
    return seqs  # 1-based ctg_id -> sequence


def parse_truth(reads_fasta):
    truth = {}
    i = 0
    with open(reads_fasta) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            i += 1
            p = line[1:].split()
            if len(p) >= 6 and p[0] == "Sim":
                truth[i] = (int(p[2]), p[3], int(p[4]), int(p[5]))
    return truth


def main():
    ref_fasta, reads_fasta = sys.argv[1], sys.argv[2]
    k = int(sys.argv[3]) if len(sys.argv) > 3 else K
    prefix = os.path.splitext(ref_fasta)[0]
    out = f"{prefix}.strictres.k{k}.txt"

    seqs = parse_ref_seqs(ref_fasta)
    truth = parse_truth(reads_fasta)

    spans = f"{prefix}.gtspans.fa"
    paf = f"{prefix}.gtspans.paf"
    n = 0
    with open(spans, "w") as fh:
        for rid, (c, strand, s, e) in truth.items():
            if c in seqs and 0 <= s < e <= len(seqs[c]) and e - s >= k:
                fh.write(f">read_{rid}\n{seqs[c][s:e]}\n")
                n += 1
    print(f"strict_resolvable: self-mapping {n} GT spans ...", file=sys.stderr)

    # --secondary=yes -p0.1 -N50: keep co-optimal secondaries so an identical
    # homolog (same score as the GT self-hit) is reported.
    subprocess.run(
        f"minimap2 -c --secondary=yes -p0.1 -N50 {ref_fasta} {spans} > {paf} 2>/dev/null",
        shell=True, check=True,
    )

    nm0 = collections.Counter()  # read_id -> # of full-length NM=0 loci
    with open(paf) as fh:
        for line in fh:
            col = line.rstrip("\n").split("\t")
            m = re.match(r"read_(\d+)$", col[0])
            if not m:
                continue
            rid = int(m.group(1))
            qlen, qs, qe = int(col[1]), int(col[2]), int(col[3])
            nm = None
            for t in col[12:]:
                if t.startswith("NM:i:"):
                    nm = int(t[5:])
            if nm == 0 and (qe - qs) >= 0.99 * qlen:
                nm0[rid] += 1

    nres = 0
    with open(out, "w") as fh:
        for rid in sorted(truth):
            res = 1 if nm0.get(rid, 0) <= 1 else 0  # only the GT at NM=0 => resolvable
            nres += res
            fh.write(f"{rid} {res}\n")
    os.remove(spans)
    os.remove(paf)
    print(f"strict_resolvable: wrote {out}  ({nres}/{len(truth)} resolvable)", file=sys.stderr)


if __name__ == "__main__":
    main()
