#!/usr/bin/env python3
"""
eval_mapping.py - Compare himapper LAshow output against ground truth
                  encoded in simulated reads FASTA headers.

FASTA header format (from simulation):
  >Sim <ref_id> <ctg_id> <strand(+/-)> <start> <end>
  (start, end are 0-based positions on the forward reference)

LAshow line format:
  >  A-id.A-ctg  B-id.B-ctg  n/c  <Ab..Ae>  x  [Bb..Be]  ~  XX.XX%  (...)

Usage:
  python eval_mapping.py <reads.fasta> <lashow_output>
  or pipe:
  LAshow reads.dam ref.dam reads.ref.las | python eval_mapping.py reads.fasta
"""

import re
import sys
from collections import Counter, defaultdict

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
OVERLAP_THRESHOLD = 0.5  # minimum reciprocal overlap fraction to call a hit correct
EDGE_MARGIN = 10000  # exclude reads overlapping within N bp of either contig end
KMER = 32  # k-mer size (matches minimaxmer K and the FastK -k32 in 0-prep.sh)


# ---------------------------------------------------------------------------
# Parse FASTA headers -> ground truth
# ---------------------------------------------------------------------------
def parse_fasta_headers(fasta_path):
    """Returns dict: read_index (1-based) -> (ctg_id, strand, start, end)"""
    truth = {}
    read_idx = 0
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if not line.startswith(">"):
                continue
            read_idx += 1
            parts = line[1:].split()
            if len(parts) < 6 or parts[0] != "Sim":
                print(
                    f"WARNING: unexpected header at read {read_idx}: {line}",
                    file=sys.stderr,
                )
                truth[read_idx] = None
                continue
            ctg_id = int(parts[2])
            strand = parts[3]
            start = int(parts[4])
            end = int(parts[5])
            truth[read_idx] = (ctg_id, strand, start, end)
    return truth


# ---------------------------------------------------------------------------
# Parse reference FASTA -> contig lengths (for end-edge exclusion)
# ---------------------------------------------------------------------------
def parse_ref_lengths(ref_fasta):
    """Returns dict: 1-based ctg_id (order of appearance) -> contig length (bp)."""
    lengths = {}
    idx = 0
    cur = 0
    with open(ref_fasta) as f:
        for line in f:
            if line.startswith(">"):
                if idx > 0:
                    lengths[idx] = cur
                idx += 1
                cur = 0
            else:
                cur += len(line.strip())
    if idx > 0:
        lengths[idx] = cur
    return lengths


# ---------------------------------------------------------------------------
# Truth-based resolvability: a read is "resolvable" iff its ground-truth origin
# span contains a genome-UNIQUE (copy-1) k-mer, i.e. a position that distinguishes
# the true locus from every homolog/repeat.  This is a property of the REFERENCE
# (independent of the tool and of read errors), unlike the diffs-of-reported-hits
# heuristic.  We read the per-position genome k-mer count straight from the FastK
# profile that minimaxmer already uses (ref.prof), via Profex.  This is the
# lightweight proxy: copy-1 presence is sufficient (not necessary) for resolvable,
# so it slightly UNDER-counts resolvable reads.
# ---------------------------------------------------------------------------
def genome_unique_positions(ref_fasta, n_ctg, K=KMER):
    """{1-based ctg_id: sorted list of positions p where ref[p:p+K] is genome-unique}.
    Returns None if the FastK profile is unavailable (truth-based report skipped)."""
    import os
    import subprocess

    prefix = os.path.splitext(ref_fasta)[0]  # ref.fasta -> ref  (FastK source root)
    if not os.path.exists(prefix + ".prof"):
        print(
            f"WARNING: FastK profile '{prefix}.prof' not found; "
            f"truth-based resolvability skipped (run FastK -k{K} -p {ref_fasta})",
            file=sys.stderr,
        )
        return None
    uniq = {}
    for c in range(1, n_ctg + 1):
        cache = f"{prefix}.uniqpos.k{K}.c{c}.txt"
        if os.path.exists(cache) and os.path.getmtime(cache) >= os.path.getmtime(
            prefix + ".prof"
        ):
            with open(cache) as fh:
                uniq[c] = [int(x) for x in fh.read().split()]
            continue
        # Profex dumps "<pos>:  <count>" per position (in order); keep count==1.
        awk = r"""$1 ~ /^[0-9]+:$/ && $2==1 {sub(/:/,"",$1); print $1}"""
        proc = subprocess.run(
            f"Profex {prefix} {c} | awk '{awk}'",
            shell=True, capture_output=True, text=True,
        )
        pos = [int(x) for x in proc.stdout.split()]
        uniq[c] = pos  # Profex emits positions in ascending order -> already sorted
        with open(cache, "w") as fh:
            fh.write("\n".join(map(str, pos)))
    return uniq


def span_has_unique(uniq_sorted, start, end, K=KMER):
    """True iff some genome-unique position p satisfies start <= p <= end-K."""
    import bisect

    hi = end - K
    if uniq_sorted is None or hi < start:
        return False
    i = bisect.bisect_left(uniq_sorted, start)
    return i < len(uniq_sorted) and uniq_sorted[i] <= hi


def load_strict_resolvable(ref_fasta, reads_fasta, K=KMER):
    """{read_id: bool} STRICT truth-based resolvability (no exact-duplicate homolog
    of the GT span).  Loads the strict_resolvable.py cache, auto-(re)generating it
    when stale.  Returns None if unavailable."""
    import os
    import subprocess

    prefix = os.path.splitext(ref_fasta)[0]
    cache = f"{prefix}.strictres.k{K}.txt"
    fresh = (
        os.path.exists(cache)
        and os.path.getmtime(cache) >= os.path.getmtime(ref_fasta)
        and os.path.getmtime(cache) >= os.path.getmtime(reads_fasta)
    )
    if not fresh:
        here = os.path.dirname(os.path.abspath(__file__))
        try:
            subprocess.run(
                f"python3 {here}/strict_resolvable.py "
                f"{ref_fasta} {reads_fasta} {K}",
                shell=True, check=True,
            )
        except Exception as exc:  # noqa: BLE001
            print(
                f"WARNING: strict resolvability precompute failed ({exc}); skipped",
                file=sys.stderr,
            )
            return None
    res = {}
    if os.path.exists(cache):
        with open(cache) as fh:
            for line in fh:
                a, b = line.split()
                res[int(a)] = b == "1"
    return res or None


# ---------------------------------------------------------------------------
# Parse LAshow output -> predictions
# ---------------------------------------------------------------------------
def parse_lashow(lashow_path_or_stream):
    """Returns dict: read_id -> list of (ctg_id, strand, ref_start, ref_end, divergence)"""
    preds = defaultdict(list)
    NUM = r"[\d,]+"
    line_re = re.compile(
        r"^[>+\- ]+"  # LAshow marks the best mapping '>' and secondary/chained ones '+';
                       # parse them all (was '^>' only, which silently dropped secondaries)
        r"(" + NUM + r")\.(\d+)\s+"
        r"(" + NUM + r")\.(\d+)\s+"
        r"([nc])\s+"
        r"[<\[]\s*(" + NUM + r")\s*\.\.\s*"
        r"(" + NUM + r")\s*[>\]]\s+x\s+"
        r"[<\[]\s*(" + NUM + r")\s*\.\.\s*"
        r"(" + NUM + r")\s*[>\]]\s+"
        r"~\s+([\d.]+)%"
        r"\s*\([\d,]+\s*x\s*[\d,]+\s*bps,\s*(" + NUM + r")\s+diffs"
    )

    def parse_num(s):
        return int(s.replace(",", ""))

    if isinstance(lashow_path_or_stream, str):
        fh = open(lashow_path_or_stream)
    else:
        fh = lashow_path_or_stream
    for line in fh:
        m = line_re.match(line)
        if not m:
            continue
        read_id = parse_num(m.group(1))
        ref_ctg = parse_num(m.group(3))
        flag = m.group(5)
        b1 = parse_num(m.group(8))
        b2 = parse_num(m.group(9))
        div = float(m.group(10))
        diffs = parse_num(m.group(11))
        strand = "+" if flag == "n" else "-"
        preds[read_id].append(
            (ref_ctg, strand, min(b1, b2), max(b1, b2), div, diffs)
        )
    if isinstance(lashow_path_or_stream, str):
        fh.close()
    return preds


# ---------------------------------------------------------------------------
# Reciprocal overlap
# ---------------------------------------------------------------------------
def reciprocal_overlap(s1, e1, s2, e2):
    inter = max(0, min(e1, e2) - max(s1, s2))
    union = max(e1, e2) - min(s1, s2)
    if union <= 0:
        return 0.0
    return inter / union


# ---------------------------------------------------------------------------
# Distinct reference loci a read maps to (same contig + overlapping = one locus)
# ---------------------------------------------------------------------------
def distinct_loci(hits):
    """Returns list of (ctg, min_diffs) for each distinct locus the read hits."""
    by_ctg = defaultdict(list)
    for ctg, strand, rs, re_, div, diffs in hits:
        by_ctg[ctg].append((rs, re_, diffs))
    loci = []
    for ctg, lst in by_ctg.items():
        lst.sort()
        ce, cd = lst[0][1], lst[0][2]
        for rs, re_, diffs in lst[1:]:
            if rs <= ce:  # overlaps the current cluster -> same locus
                ce = max(ce, re_)
                cd = min(cd, diffs)
            else:
                loci.append((ctg, cd))
                ce, cd = re_, diffs
        loci.append((ctg, cd))
    return loci


# ---------------------------------------------------------------------------
# Comprehensive evaluate
# ---------------------------------------------------------------------------
def evaluate(
    truth, preds, ctg_len=None, gt_uniq=None, gt_strict=None,
    overlap_thr=OVERLAP_THRESHOLD,
):
    pct = lambda n, d: f"{100*n/d:.1f}%" if d > 0 else "N/A"
    ratio = lambda n, d: f"{n/d:.4f}" if d > 0 else "N/A"
    skipped_edge = 0
    total = 0
    unmapped_ids = []
    wrong_ids = []

    # per-read best-hit
    mapped = 0
    fully_correct = 0

    # per-mapping
    total_mappings = 0
    m_fc = 0
    m_fp = 0
    reads_w_correct = 0

    # multiplicity & contig pattern
    n_maps_dist = Counter()
    ctg_sets = Counter()

    # wrong breakdown
    wrong_detail = Counter()
    wrong_detail_ids = defaultdict(list)

    # resolvability: identical-sequence ambiguity (excused) vs >=1 base difference
    resolvable = 0
    res_correct = 0
    unresolvable = 0
    unres_correct = 0
    res_wrong = Counter()
    res_wrong_ids = defaultdict(list)
    # truth-based resolvability (genome-unique k-mer in the GT span; tool-independent)
    t_res = 0
    t_res_correct = 0
    t_res_wrong_ids = []
    # strict truth-based resolvability (no exact-duplicate homolog of the GT span)
    s_res = 0
    s_res_correct = 0
    s_res_wrong_ids = []

    for read_id in sorted(truth.keys()):
        gt = truth[read_id]
        if gt is None:
            continue
        gt_ctg, gt_strand, gt_start, gt_end = gt
        clen = ctg_len.get(gt_ctg) if ctg_len else None
        near_start = gt_start < EDGE_MARGIN
        near_end = clen is not None and gt_end > clen - EDGE_MARGIN
        crosses_boundary = gt_end < gt_start  # wrap read across a contig seam
        if near_start or near_end or crosses_boundary:
            skipped_edge += 1
            continue
        total += 1
        hits = preds.get(read_id, [])
        n_maps_dist[len(hits)] += 1

        if not hits:
            unmapped_ids.append(read_id)
            ctg_sets[frozenset()] += 1
            continue

        hit_ctgs = frozenset(h[0] for h in hits)
        ctg_sets[hit_ctgs] += 1
        mapped += 1

        # -- per-read best-hit --
        def score(h):
            ctg, strand, rs, re, div, diffs = h
            ov = reciprocal_overlap(gt_start, gt_end, rs, re) if ctg == gt_ctg else 0.0
            return (strand == gt_strand, ov, -div)

        best = max(hits, key=score)
        b_ctg, b_strand, b_start, b_end, b_div, b_diffs = best
        ov = (
            reciprocal_overlap(gt_start, gt_end, b_start, b_end)
            if b_ctg == gt_ctg
            else 0.0
        )
        if b_strand == gt_strand and ov >= overlap_thr:
            fully_correct += 1
        else:
            wrong_ids.append(read_id)
            if b_ctg == gt_ctg:
                category = (
                    "strand_ok, pos_wrong (same ctg)"
                    if b_strand == gt_strand
                    else "strand_wrong (same ctg)"
                )
            else:
                category = (
                    "strand_ok, wrong ctg"
                    if b_strand == gt_strand
                    else "strand_wrong, wrong ctg"
                )
            wrong_detail[category] += 1
            wrong_detail_ids[category].append(read_id)

        # -- resolvability: identical-sequence tie vs >=1 base difference --
        loci = distinct_loci(hits)
        d_min = min(d for _, d in loci)
        n_min = sum(1 for _, d in loci if d == d_min)
        gt_diffs = min(
            (
                h[5]
                for h in hits
                if h[0] == gt_ctg
                and reciprocal_overlap(gt_start, gt_end, h[2], h[3]) >= 0.3
            ),
            default=None,
        )
        is_correct = b_strand == gt_strand and ov >= overlap_thr

        # truth-based resolvability: is the GT origin span distinguishable at all?
        if gt_uniq is not None and span_has_unique(
            gt_uniq.get(gt_ctg), gt_start, gt_end
        ):
            t_res += 1
            if is_correct:
                t_res_correct += 1
            else:
                t_res_wrong_ids.append(read_id)

        if gt_strict is not None and gt_strict.get(read_id):
            s_res += 1
            if is_correct:
                s_res_correct += 1
            else:
                s_res_wrong_ids.append(read_id)

        if n_min >= 2:  # >=2 distinct loci tie at min diffs => identical sequence
            unresolvable += 1
            if is_correct:
                unres_correct += 1
        else:  # unique min-diff locus => >=1 base distinguishes the true origin
            resolvable += 1
            if is_correct:
                res_correct += 1
            elif gt_diffs is None:
                res_wrong["GT locus not found"] += 1
                res_wrong_ids["GT locus not found"].append(read_id)
            elif gt_diffs <= d_min:
                key = "GT best by diffs but not chosen (selection)"
                res_wrong[key] += 1
                res_wrong_ids[key].append(read_id)
            else:
                res_wrong["read fits other locus better"] += 1
                res_wrong_ids["read fits other locus better"].append(read_id)

        # -- per-mapping --
        any_correct = False
        for h in hits:
            ctg, strand, rs, re, div, diffs = h
            total_mappings += 1
            ov_m = (
                reciprocal_overlap(gt_start, gt_end, rs, re) if ctg == gt_ctg else 0.0
            )
            if strand == gt_strand and ov_m >= overlap_thr:
                m_fc += 1
                any_correct = True
            else:
                m_fp += 1
        if any_correct:
            reads_w_correct += 1

    mapping_precision = m_fc / total_mappings if total_mappings > 0 else None
    mapping_recall = reads_w_correct / total if total > 0 else None
    mapping_f1 = None
    if mapping_precision is not None and mapping_recall is not None:
        denom = mapping_precision + mapping_recall
        mapping_f1 = (
            (2 * mapping_precision * mapping_recall / denom) if denom > 0 else 0.0
        )

    multimapped_reads = sum(cnt for k, cnt in n_maps_dist.items() if k > 1)
    cross_contig_reads = sum(cnt for cs, cnt in ctg_sets.items() if len(cs) > 1)

    # --- output ---
    print("=== Mapping Accuracy ===")
    print(f"  Number of reads     : {total}")
    edge_where = "either end" if ctg_len else "start"
    print(
        f"  Edge-excluded reads : {skipped_edge}  (within {EDGE_MARGIN}bp of contig {edge_where})"
    )
    print()

    print("  [Read-level best hit]")
    print(f"    Mapped            : {mapped}  ({pct(mapped, total)})")
    print(
        f"    Unmapped          : {len(unmapped_ids)}  ({pct(len(unmapped_ids), total)})"
    )
    print(f"    Correct hit       : {fully_correct}  ({pct(fully_correct, total)})")
    print(f"    Wrong hit         : {len(wrong_ids)}  ({pct(len(wrong_ids), total)})")
    import os as _os
    if _os.getenv("DUMP_WRONG"):
        with open(_os.getenv("DUMP_WRONG"), "w") as _f:
            _f.write(" ".join(map(str, sorted(wrong_ids))))
    print(f"      (same contig, same strand, overlap >= {overlap_thr:.0%})")
    print()

    print("  [Mapping-level]")
    print(f"    Predicted         : {total_mappings}")
    if total_mappings > 0:
        print(f"    True positive     : {m_fc}  ({pct(m_fc, total_mappings)})")
        print(f"    False positive    : {m_fp}  ({pct(m_fp, total_mappings)})")
    print(f"    Recall            : {ratio(reads_w_correct, total)}")
    print(f"    Precision         : {ratio(m_fc, total_mappings)}")
    if mapping_f1 is not None:
        print(f"    F1                : {mapping_f1:.4f}")
    print()

    print("  [Ambiguity diagnostics]")
    print(
        f"    Multi-mapped reads: {multimapped_reads}  ({pct(multimapped_reads, total)})"
    )
    print(
        f"    Cross-contig reads: {cross_contig_reads}  ({pct(cross_contig_reads, total)})"
    )
    print()

    print("  [Resolvability]  (identical-sequence ambiguity excused; needs all mappings, run WITHOUT -1)")
    print(
        f"    Unresolvable ties : {unresolvable}  ({pct(unresolvable, total)})"
        f"   [>=2 loci equal at min diffs = identical sequence]"
    )
    print(f"      best=GT (tie-break): {unres_correct}")
    print(
        f"    Resolvable        : {resolvable}  ({pct(resolvable, total)})"
        f"   [unique min-diff locus = >=1 base distinguishes truth]"
    )
    print(f"      Correct          : {res_correct}  ({pct(res_correct, resolvable)})")
    print(f"      Wrong (fixable)  : {resolvable - res_correct}")
    for cat, cnt in res_wrong.most_common():
        print(f"        {cat:42s}: {cnt}")
        print(f"          first IDs: {res_wrong_ids[cat][:10]}")
    print(
        f"    >> Resolvable accuracy: {ratio(res_correct, resolvable)}"
        f"   (excludes {unresolvable} unresolvable ties)"
    )
    print()

    if gt_uniq is not None:
        print(
            "  [Resolvability - truth-based]"
            "  (GT span has a genome-unique k-mer; tool-independent)"
        )
        print(
            f"    Resolvable (truth): {t_res}  ({pct(t_res, total)})"
            f"   [>=1 copy-1 k-mer in the true origin span]"
        )
        print(f"      Correct          : {t_res_correct}  ({pct(t_res_correct, t_res)})")
        print(f"      Wrong            : {t_res - t_res_correct}")
        print(f"        first IDs: {sorted(t_res_wrong_ids)[:15]}")
        print(
            f"    >> Resolvable accuracy (truth): {ratio(t_res_correct, t_res)}"
            f"   (lightweight copy-1 proxy; conservative)"
        )
        print()

    if gt_strict is not None:
        print(
            "  [Resolvability - strict]"
            "  (GT span has no exact-duplicate homolog; tool-independent)"
        )
        print(
            f"    Resolvable (strict): {s_res}  ({pct(s_res, total)})"
            f"   [no 0-diff full-length homolog of the true origin span]"
        )
        print(f"      Correct          : {s_res_correct}  ({pct(s_res_correct, s_res)})")
        print(f"      Wrong            : {s_res - s_res_correct}")
        print(f"        first IDs: {sorted(s_res_wrong_ids)[:15]}")
        print(
            f"    >> Resolvable accuracy (strict): {ratio(s_res_correct, s_res)}"
            f"   (rigorous: self-aligned reference)"
        )
        print()

    if wrong_detail:
        print(f"  [Wrong reads breakdown ({len(wrong_ids)} reads)]")
        for cat, cnt in sorted(wrong_detail.items(), key=lambda x: -x[1]):
            print(f"    {cat:40s} : {cnt}")
            print(f"      first read IDs: {wrong_detail_ids[cat][:10]}")
        print()

    if unmapped_ids[:10]:
        print(f"  First unmapped read IDs: {unmapped_ids[:10]}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Evaluate himapper LAshow output against ground truth."
    )
    parser.add_argument("reads_fasta", help="Reads FASTA with ground-truth headers")
    parser.add_argument(
        "lashow", nargs="?", default=None, help="LAshow output file (default: stdin)"
    )
    parser.add_argument(
        "-r",
        "--ref",
        default=None,
        help="Reference FASTA; if given, also exclude reads within "
        "EDGE_MARGIN bp of the contig END (not only the start)",
    )
    args = parser.parse_args()

    fasta_path = args.reads_fasta
    lashow_source = args.lashow if args.lashow else sys.stdin
    truth = parse_fasta_headers(fasta_path)
    preds = parse_lashow(lashow_source)
    ctg_len = parse_ref_lengths(args.ref) if args.ref else None
    gt_uniq = (
        genome_unique_positions(args.ref, len(ctg_len)) if args.ref and ctg_len else None
    )
    gt_strict = load_strict_resolvable(args.ref, fasta_path) if args.ref else None
    evaluate(truth, preds, ctg_len, gt_uniq, gt_strict)
