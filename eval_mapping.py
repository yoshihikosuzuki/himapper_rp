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

import sys
import re
import statistics
from collections import defaultdict, Counter

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
OVERLAP_THRESHOLD = 0.5   # minimum reciprocal overlap fraction to call a hit correct
EDGE_MARGIN       = 10000 # exclude reads that overlap the first N bp of a contig


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
# Parse LAshow output -> predictions
# ---------------------------------------------------------------------------
def parse_lashow(lashow_path_or_stream):
    """Returns dict: read_id -> list of (ctg_id, strand, ref_start, ref_end, divergence)"""
    preds = defaultdict(list)
    NUM  = r'[\d,]+'
    line_re = re.compile(
        r'^>\s+'
        r'(' + NUM + r')\.(\d+)\s+'
        r'(' + NUM + r')\.(\d+)\s+'
        r'([nc])\s+'
        r'[<\[]\s*(' + NUM + r')\s*\.\.\s*'
        r'(' + NUM + r')\s*[>\]]\s+x\s+'
        r'[<\[]\s*(' + NUM + r')\s*\.\.\s*'
        r'(' + NUM + r')\s*[>\]]\s+'
        r'~\s+([\d.]+)%'
    )
    def parse_num(s): return int(s.replace(',', ''))

    if isinstance(lashow_path_or_stream, str):
        fh = open(lashow_path_or_stream)
    else:
        fh = lashow_path_or_stream
    for line in fh:
        m = line_re.match(line)
        if not m: continue
        read_id  = parse_num(m.group(1))
        ref_ctg  = parse_num(m.group(3))
        flag     = m.group(5)
        b1       = parse_num(m.group(8))
        b2       = parse_num(m.group(9))
        div      = float(m.group(10))
        strand   = '+' if flag == 'n' else '-'
        preds[read_id].append((ref_ctg, strand, min(b1,b2), max(b1,b2), div))
    if isinstance(lashow_path_or_stream, str):
        fh.close()
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
            ctg, strand, rs, re, div = h
            ov = reciprocal_overlap(gt_start, gt_end, rs, re) if ctg == gt_ctg else 0.0
            return (strand == gt_strand, ov, -div)

        best = max(hits, key=score)
        b_ctg, b_strand, b_start, b_end, b_div = best
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
            ctg, strand, rs, re, div = h
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
        description='Evaluate himapper LAshow output against ground truth.')
    parser.add_argument('reads_fasta', help='Reads FASTA with ground-truth headers')
    parser.add_argument('lashow', nargs='?', default=None,
                        help='LAshow output file (default: stdin)')
    args = parser.parse_args()

    fasta_path = args.reads_fasta
    lashow_source = args.lashow if args.lashow else sys.stdin
    truth = parse_fasta_headers(fasta_path)
    preds = parse_lashow(lashow_source)
    evaluate(truth, preds)
