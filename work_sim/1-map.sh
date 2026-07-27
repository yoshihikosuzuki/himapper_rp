#!/bin/bash
# ml daligner

echo "=== minimaxmer ==="
{ time {
  ../minimaxmer -T1 ref.dam ref.prof ref.reads.prof .ref.minimax.kpos
  ../minimaxmer -T1 reads.dam reads.ref.prof reads.prof .reads.minimax.kpos
  ../himapper -v -T1 \
    ref.dam .ref.minimax.kpos \
    reads.dam .reads.minimax.kpos
}; } 2>&1
mv reads.ref.las reads.ref.minimax.las
LAshow reads.dam ref.dam reads.ref.minimax.las | python3 scripts/eval_mapping.py reads.fasta --ref ref.fasta >mapping_results.minimax.txt
cat mapping_results.minimax.txt


echo "=== minimap2 ==="
# minimap2 -x map-hifi ref.fasta reads.fasta >minimap2.paf
# minimap2 -x map-hifi --secondary=no ref.fasta reads.fasta >minimap2_nosec.paf
python3 scripts/eval_minimap2.py ref.fasta reads.fasta minimap2_nosec.paf >mapping_results.mm2.txt
cat mapping_results.mm2.txt

echo "=== winnowmap2 ==="
# winnowmap -x map-pb ref.fasta reads.fasta >winnowmap2.paf
# winnowmap -x map-pb --secondary=no ref.fasta reads.fasta >winnowmap2_nosec.paf
python3 scripts/eval_minimap2.py ref.fasta reads.fasta winnowmap_nosec.paf >mapping_results.wm2.txt
cat mapping_results.wm2.txt
