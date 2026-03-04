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
from collections import defaultdict

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
OVERLAP_THRESHOLD = 0.5   # minimum reciprocal overlap fraction to call a hit correct
POS_TOLERANCE     = 500   # bp tolerance for start/end positions (informational)
EDGE_MARGIN       = 10000 # exclude reads that overlap the first N bp of a contig


# ---------------------------------------------------------------------------
# Parse FASTA headers -> ground truth
# ---------------------------------------------------------------------------
def parse_fasta_headers(fasta_path):
    """
    Returns dict: read_index (1-based) -> (ctg_id, strand, start, end)
    """
    truth = {}
    read_idx = 0
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if not line.startswith('>'):
                continue
            read_idx += 1
            # >Sim ref_id ctg_id strand start end
            parts = line[1:].split()
            if len(parts) < 6 or parts[0] != 'Sim':
                print(f"WARNING: unexpected header at read {read_idx}: {line}", file=sys.stderr)
                truth[read_idx] = None
                continue
            ctg_id = int(parts[2])
            strand = parts[3]          # '+' or '-'
            start  = int(parts[4])
            end    = int(parts[5])
            truth[read_idx] = (ctg_id, strand, start, end)
    return truth


# ---------------------------------------------------------------------------
# Parse LAshow output -> predictions
# ---------------------------------------------------------------------------
def parse_lashow(lashow_path_or_stream):
    """
    Returns dict: read_id (1-based int) -> list of (ctg_id, strand, ref_start, ref_end, divergence)
    Multiple mappings per read are kept; caller can choose best.
    """
    preds = defaultdict(list)

    # Regex for a LAshow data line, handling comma-formatted numbers and
    # variable whitespace (LAshow pads numbers in fixed-width columns so
    # there can be spaces both before and after "..").
    # Note: long lines may wrap after "~XX.XX%  (" — we only need up to %.
    NUM  = r'[\d,]+'
    line_re = re.compile(
        r'^>\s+'
        r'(' + NUM + r')\.(\d+)\s+'         # A-id, A-ctg
        r'(' + NUM + r')\.(\d+)\s+'         # B-id, B-ctg
        r'([nc])\s+'                         # strand flag n/c
        r'[<\[]\s*(' + NUM + r')\s*\.\.\s*' # A-start (skip A coords)
        r'(' + NUM + r')\s*[>\]]\s+x\s+'    # A-end
        r'[<\[]\s*(' + NUM + r')\s*\.\.\s*' # B-start
        r'(' + NUM + r')\s*[>\]]\s+'        # B-end
        r'~\s+([\d.]+)%'                    # divergence
    )

    def parse_num(s):
        return int(s.replace(',', ''))

    if isinstance(lashow_path_or_stream, str):
        fh = open(lashow_path_or_stream)
    else:
        fh = lashow_path_or_stream

    for line in fh:
        m = line_re.match(line)
        if not m:
            continue
        read_id  = parse_num(m.group(1))
        ref_ctg  = parse_num(m.group(3))
        flag     = m.group(5)           # 'n' or 'c'
        b1       = parse_num(m.group(8))
        b2       = parse_num(m.group(9))
        div      = float(m.group(10))

        strand       = '+' if flag == 'n' else '-'
        ref_start    = min(b1, b2)   # forward coords in both n and c cases
        ref_end      = max(b1, b2)

        preds[read_id].append((ref_ctg, strand, ref_start, ref_end, div))

    if isinstance(lashow_path_or_stream, str):
        fh.close()

    return preds


# ---------------------------------------------------------------------------
# Reciprocal overlap
# ---------------------------------------------------------------------------
def reciprocal_overlap(s1, e1, s2, e2):
    """Fraction of intersection over union of two intervals."""
    inter = max(0, min(e1, e2) - max(s1, s2))
    union = max(e1, e2) - min(s1, s2)
    if union <= 0:
        return 0.0
    return inter / union


# ---------------------------------------------------------------------------
# Evaluate
# ---------------------------------------------------------------------------
def evaluate(truth, preds, overlap_thr=OVERLAP_THRESHOLD):
    total         = len(truth)
    mapped        = 0          # has at least one mapping
    strand_correct = 0         # best hit has correct strand
    pos_correct    = 0         # best hit has overlap >= threshold (correct strand implied)
    fully_correct  = 0         # correct strand AND position

    start_errs = []
    end_errs   = []
    unmapped   = []
    skipped_edge = 0           # reads within EDGE_MARGIN of contig start

    for read_id in sorted(truth.keys()):
        gt = truth[read_id]
        if gt is None:
            total -= 1
            continue
        gt_ctg, gt_strand, gt_start, gt_end = gt

        # Exclude reads that overlap the first EDGE_MARGIN bp of the contig.
        if gt_start < EDGE_MARGIN:
            total -= 1
            skipped_edge += 1
            continue

        hits = preds.get(read_id, [])
        if not hits:
            unmapped.append(read_id)
            continue
        mapped += 1

        # choose best hit: prefer correct strand, then highest overlap, then lowest divergence
        def score(h):
            ctg, strand, rs, re, div = h
            ov = reciprocal_overlap(gt_start, gt_end, rs, re) if ctg == gt_ctg else 0.0
            return (strand == gt_strand, ov, -div)

        best = max(hits, key=score)
        b_ctg, b_strand, b_start, b_end, b_div = best

        if b_strand == gt_strand:
            strand_correct += 1

        ov = reciprocal_overlap(gt_start, gt_end, b_start, b_end) if b_ctg == gt_ctg else 0.0
        if ov >= overlap_thr:
            pos_correct += 1

        if b_strand == gt_strand and ov >= overlap_thr:
            fully_correct += 1
            start_errs.append(abs(b_start - gt_start))
            end_errs.append(abs(b_end - gt_end))

    print(f"=== Mapping Accuracy ===")
    print(f"  Edge-excluded reads : {skipped_edge}  (start < {EDGE_MARGIN} bp from contig start)")
    print(f"  Total reads         : {total}")
    print(f"  Mapped              : {mapped}  ({100*mapped/total:.1f}%)")
    print(f"  Unmapped            : {len(unmapped)}  ({100*len(unmapped)/total:.1f}%)")
    print(f"  Strand correct      : {strand_correct}  ({100*strand_correct/total:.1f}%)")
    print(f"  Position correct    : {pos_correct}  ({100*pos_correct/total:.1f}%)")
    print(f"    (overlap_thr = {overlap_thr:.0%}, same contig)")
    print(f"  Fully correct       : {fully_correct}  ({100*fully_correct/total:.1f}%)")
    print()
    if start_errs:
        import statistics
        print(f"  Among fully correct ({fully_correct} reads):")
        print(f"    Start err  median={statistics.median(start_errs):6.0f}  mean={sum(start_errs)/len(start_errs):6.0f}  max={max(start_errs)}")
        print(f"    End   err  median={statistics.median(end_errs):6.0f}  mean={sum(end_errs)/len(end_errs):6.0f}  max={max(end_errs)}")

    if unmapped[:10]:
        print(f"\n  First unmapped read IDs: {unmapped[:10]}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    fasta_path = sys.argv[1]

    if len(sys.argv) >= 3:
        lashow_source = sys.argv[2]
    else:
        lashow_source = sys.stdin

    truth = parse_fasta_headers(fasta_path)
    preds = parse_lashow(lashow_source)

    evaluate(truth, preds)
