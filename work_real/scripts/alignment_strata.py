#!/usr/bin/env python3
"""Estimate mapping difficulty by self-mapping mapping-truth reference spans.

For each read, extract its assigned reference span and map it back to the diploid
panel while retaining secondary alignments.

Output:
  <read_name> <level> <n_nm0> <nearest_other_NM>

mapping_truth.tsv columns: read_name, hap, contig, start, end

Usage:
  alignment_strata.py <panel.fasta> <mapping_truth.tsv> <out.tsv>
      [--threads 8] [--near 2] [--med 10]
"""
import argparse
import collections
import os
import subprocess
import sys


def parse_ref(fa):
    """{contig_name: sequence}.  Name = first whitespace token of the header."""
    seqs, name, cur = {}, None, []
    with open(fa) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(cur)
                name = line[1:].split()[0]
                cur = []
            else:
                cur.append(line.strip())
        if name is not None:
            seqs[name] = "".join(cur)
    return seqs


def parse_mapping_truth(mapping_truth_tsv):
    """{read_name: (ctg, start, end)}."""
    mapping_truth = {}
    with open(mapping_truth_tsv) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                mapping_truth[p[0]] = (p[2], int(p[3]), int(p[4]))
    return mapping_truth


def selfmap_full_length_nm(panel, seqs, mapping_truth, out, threads):
    """Cut assigned spans, self-map to the panel, and return full-length-hit NMs."""
    spans, paf = out + ".spans.fa", out + ".spans.paf"
    id2name = {}
    with open(spans, "w") as fh:
        for i, (name, (c, s, e)) in enumerate(mapping_truth.items()):
            if c in seqs and 0 <= s < e <= len(seqs[c]):
                sid = f"span_{i}"
                id2name[sid] = name
                fh.write(f">{sid}\n{seqs[c][s:e]}\n")
    print(
        f"alignment_strata: self-mapping {len(id2name)} assigned spans ...",
        file=sys.stderr)
    subprocess.run(
        f"minimap2 -t {threads} -c --secondary=yes -p0.1 -N50 {panel} {spans} > {paf} 2>/dev/null",
        shell=True, check=True)

    nms = collections.defaultdict(list)  # read_name -> [NM of full-length hits]
    with open(paf) as fh:
        for line in fh:
            col = line.rstrip("\n").split("\t")
            sid = col[0]
            if sid not in id2name:
                continue
            qlen, qs, qe = int(col[1]), int(col[2]), int(col[3])
            nm = next((int(t[5:]) for t in col[12:] if t.startswith("NM:i:")), None)
            if nm is not None and (qe - qs) >= 0.99 * qlen:
                nms[id2name[sid]].append(nm)
    os.remove(spans)
    os.remove(paf)
    return {name: sorted(v) for name, v in nms.items()}, id2name


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("panel")
    ap.add_argument("mapping_truth")
    ap.add_argument("out")
    ap.add_argument("--threads", type=int, default=8,
                    help="number of minimap2 threads [8]")
    ap.add_argument("--near", type=int, default=2, help="hard NM threshold")
    ap.add_argument("--med", type=int, default=10, help="med NM threshold")
    args = ap.parse_args()

    seqs = parse_ref(args.panel)
    mapping_truth = parse_mapping_truth(args.mapping_truth)
    nms, id2name = selfmap_full_length_nm(
        args.panel, seqs, mapping_truth, args.out, args.threads)

    counts = collections.Counter()
    with open(args.out, "w") as fh:
        for _sid, name in id2name.items():
            v = nms.get(name, [])
            n0 = sum(1 for x in v if x == 0)
            others = v if n0 == 0 else v[1:]  # drop one NM=0 self when present
            nearest = min(others) if others else 999999
            if n0 == 0:
                lvl = "no_map"
            elif n0 >= 2:
                lvl = "dup"
            elif nearest <= args.near:
                lvl = "hard"
            elif nearest <= args.med:
                lvl = "med"
            else:
                lvl = "easy"
            counts[lvl] += 1
            fh.write(f"{name}\t{lvl}\t{n0}\t{nearest if nearest < 999999 else -1}\n")
    print("alignment_strata: " + "  ".join(
        f"{k}={counts[k]}"
        for k in ["easy", "med", "hard", "dup", "no_map"]),
        file=sys.stderr)


if __name__ == "__main__":
    main()
