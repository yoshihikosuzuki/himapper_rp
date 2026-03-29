#!/usr/bin/env python3
"""
eval_minimap2.py - Compare minimap2/winnowmap PAF output against ground truth
                   encoded in simulated reads FASTA headers.

FASTA header format (from simulation):
  >Sim <ref_id> <ctg_id> <strand(+/-)> <start> <end>
  (start, end are 0-based positions on the forward reference)

Reads must be renamed so that the query name in PAF is "read_N" (1-based).
For example:  work/reads.renamed.fasta

Usage:
  python eval_minimap2.py ref.fasta reads.fasta alignment.paf
"""

import sys
import re
import statistics
from collections import defaultdict, Counter

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
OVERLAP_THRESHOLD = 0.5
EDGE_MARGIN       = 10000


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
            if not line.startswith('>'):
                continue
            read_idx += 1
            parts = line[1:].split()
            if len(parts) < 6 or parts[0] != 'Sim':
                print(f"WARNING: unexpected header at read {read_idx}: {line}", file=sys.stderr)
                truth[read_idx] = None
                continue
            ctg_id = int(parts[2])
            strand = parts[3]
            start  = int(parts[4])
            end    = int(parts[5])
            truth[read_idx] = (ctg_id, strand, start, end)
    return truth


# ---------------------------------------------------------------------------
# Parse reference FASTA -> contig name-to-id mapping
# ---------------------------------------------------------------------------
def parse_ref_names(ref_fasta):
    """Returns dict: contig_name -> 1-based ctg_id (order of appearance)."""
    name_to_id = {}
    idx = 0
    with open(ref_fasta) as f:
        for line in f:
            if line.startswith('>'):
                idx += 1
                name = line[1:].split()[0]
                name_to_id[name] = idx
    return name_to_id


# ---------------------------------------------------------------------------
# Parse PAF -> predictions
# ---------------------------------------------------------------------------
def parse_paf(paf_path, name_to_id):
    """Returns dict: read_id (int) -> list of (ctg_id, strand, ref_start, ref_end)"""
    preds = defaultdict(list)
    with open(paf_path) as f:
        for line in f:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 12: continue
            qname  = cols[0]
            strand = cols[4]
            tname  = cols[5]
            tstart = int(cols[7])
            tend   = int(cols[8])
            # extract read id
            m = re.match(r'read_(\d+)', qname)
            if not m:
                continue
            read_id = int(m.group(1))
            ctg_id = name_to_id.get(tname)
            if ctg_id is None:
                print(f"WARNING: unknown target {tname}", file=sys.stderr)
                continue
            preds[read_id].append((ctg_id, strand, tstart, tend))
    return preds


# ---------------------------------------------------------------------------
# Reciprocal overlap
# ---------------------------------------------------------------------------
def reciprocal_overlap(s1, e1, s2, e2):
    inter = max(0, min(e1, e2) - max(s1, s2))
    union = max(e1, e2) - min(s1, s2)
    if union <= 0: return 0.0
    return inter / union


# ---------------------------------------------------------------------------
# Comprehensive evaluate
# ---------------------------------------------------------------------------
def evaluate(truth, preds, overlap_thr=OVERLAP_THRESHOLD):
    pct = lambda n, d: f"{100*n/d:.1f}%" if d > 0 else "N/A"
    skipped_edge = 0; total = 0
    unmapped_ids = []; wrong_ids = []

    # per-read best-hit
    mapped = 0; strand_correct = 0; pos_correct = 0; fully_correct = 0
    start_errs = []; end_errs = []

    # per-mapping
    total_mappings = 0; m_fc = 0; reads_w_correct = 0
    pm_start_errs = []; pm_end_errs = []

    # multiplicity & contig pattern
    n_maps_dist = Counter()
    ctg_sets    = Counter()

    # wrong breakdown
    wrong_detail = Counter()

    for read_id in sorted(truth.keys()):
        gt = truth[read_id]
        if gt is None: continue
        gt_ctg, gt_strand, gt_start, gt_end = gt
        if gt_start < EDGE_MARGIN:
            skipped_edge += 1; continue
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
            ctg, strand, rs, re = h
            ov = reciprocal_overlap(gt_start, gt_end, rs, re) if ctg == gt_ctg else 0.0
            return (strand == gt_strand, ov)

        best = max(hits, key=score)
        b_ctg, b_strand, b_start, b_end = best
        if b_strand == gt_strand: strand_correct += 1
        ov = reciprocal_overlap(gt_start, gt_end, b_start, b_end) if b_ctg == gt_ctg else 0.0
        if ov >= overlap_thr: pos_correct += 1
        if b_strand == gt_strand and ov >= overlap_thr:
            fully_correct += 1
            start_errs.append(abs(b_start - gt_start))
            end_errs.append(abs(b_end - gt_end))
        else:
            wrong_ids.append(read_id)
            if b_ctg == gt_ctg:
                wrong_detail['strand_ok, pos_wrong (same ctg)' if b_strand == gt_strand
                             else 'strand_wrong (same ctg)'] += 1
            else:
                wrong_detail['strand_ok, wrong ctg' if b_strand == gt_strand
                             else 'strand_wrong, wrong ctg'] += 1

        # -- per-mapping --
        any_correct = False
        for h in hits:
            ctg, strand, rs, re = h
            total_mappings += 1
            ov_m = reciprocal_overlap(gt_start, gt_end, rs, re) if ctg == gt_ctg else 0.0
            if strand == gt_strand and ov_m >= overlap_thr:
                m_fc += 1; any_correct = True
                pm_start_errs.append(abs(rs - gt_start))
                pm_end_errs.append(abs(re - gt_end))
        if any_correct: reads_w_correct += 1

    # --- output ---
    print(f"=== Mapping Accuracy ===")
    print(f"  Edge-excluded reads : {skipped_edge}  (gt_start < {EDGE_MARGIN} bp)")
    print()

    print(f"  [Per-read, best hit]")
    print(f"    Total reads       : {total}")
    print(f"    Mapped            : {mapped}  ({pct(mapped, total)})")
    print(f"    Unmapped          : {len(unmapped_ids)}  ({pct(len(unmapped_ids), total)})")
    print(f"    Strand correct    : {strand_correct}  ({pct(strand_correct, total)})")
    print(f"    Position correct  : {pos_correct}  ({pct(pos_correct, total)})")
    print(f"      (overlap_thr = {overlap_thr:.0%}, same contig)")
    print(f"    Fully correct     : {fully_correct}  ({pct(fully_correct, total)})")
    print(f"    Wrong (mapped,!FC): {len(wrong_ids)}  ({pct(len(wrong_ids), total)})")
    if start_errs:
        print(f"    Start err (FC)    : median={statistics.median(start_errs):6.0f}  "
              f"mean={sum(start_errs)/len(start_errs):6.0f}  max={max(start_errs)}")
        print(f"    End   err (FC)    : median={statistics.median(end_errs):6.0f}  "
              f"mean={sum(end_errs)/len(end_errs):6.0f}  max={max(end_errs)}")
    print()

    print(f"  [Per-mapping]")
    print(f"    Total mappings           : {total_mappings}")
    if total_mappings > 0:
        print(f"    Fully correct mappings   : {m_fc}  ({pct(m_fc, total_mappings)} of mappings)")
    print(f"    Reads w/ >=1 correct hit : {reads_w_correct}  ({pct(reads_w_correct, total)} of reads)")
    if pm_start_errs:
        print(f"    Start err (FC mappings)  : median={statistics.median(pm_start_errs):6.0f}  "
              f"mean={sum(pm_start_errs)/len(pm_start_errs):6.0f}  max={max(pm_start_errs)}")
        print(f"    End   err (FC mappings)  : median={statistics.median(pm_end_errs):6.0f}  "
              f"mean={sum(pm_end_errs)/len(pm_end_errs):6.0f}  max={max(pm_end_errs)}")
    print()

    print(f"  [Mapping multiplicity]")
    for k in sorted(n_maps_dist.keys()):
        print(f"    {k} mapping(s) : {n_maps_dist[k]:5d} reads")
    print()

    print(f"  [Contig hit pattern]")
    for cs in sorted(ctg_sets.keys(), key=lambda s: (len(s), sorted(s) if s else [])):
        label = "unmapped" if not cs else "+".join(f"ctg{c}" for c in sorted(cs))
        print(f"    {label:20s} : {ctg_sets[cs]:5d}  ({pct(ctg_sets[cs], total)})")
    print()

    if wrong_detail:
        print(f"  [Wrong reads breakdown ({len(wrong_ids)} reads)]")
        for cat, cnt in sorted(wrong_detail.items(), key=lambda x: -x[1]):
            print(f"    {cat:40s} : {cnt}")
        print()

    if unmapped_ids[:10]:
        print(f"  First unmapped read IDs: {unmapped_ids[:10]}")
    if wrong_ids[:10]:
        print(f"  First wrong    read IDs: {wrong_ids[:10]}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Evaluate minimap2/winnowmap PAF output against ground truth.')
    parser.add_argument('ref_fasta', help='Reference FASTA (for contig name -> id mapping)')
    parser.add_argument('reads_fasta', help='Reads FASTA with ground-truth headers')
    parser.add_argument('paf', help='PAF alignment file')
    args = parser.parse_args()

    name_to_id = parse_ref_names(args.ref_fasta)
    truth  = parse_fasta_headers(args.reads_fasta)
    preds  = parse_paf(args.paf, name_to_id)
    evaluate(truth, preds)
