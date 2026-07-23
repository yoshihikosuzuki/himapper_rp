#!/bin/bash
# ml daligner

echo "=== modimizer ==="
{ time {
  ../modimizer ref.dam ref.mod.kpos
  ../himapper -v -T1 ref.dam ref.mod.kpos reads.dam all
}; } 2>&1
mv reads.ref.las reads.ref.mod.las
LAshow reads.dam ref.dam reads.ref.mod.las | python3 scripts/eval_mapping.py reads.fasta --ref ref.fasta >mapping_results.mod.txt
cat mapping_results.mod.txt

echo "=== copy_minimizer ==="
{ time {
  ../copy_minimizer ref.dam ref.prof ref.min.kpos
  ../himapper -v -T1 ref.dam ref.min.kpos reads.dam all
}; } 2>&1
mv reads.ref.las reads.ref.min.las
LAshow reads.dam ref.dam reads.ref.min.las | python3 scripts/eval_mapping.py reads.fasta --ref ref.fasta >mapping_results.min.txt
cat mapping_results.min.txt

echo "=== minimaxmer ==="
{ time {
  ../minimaxmer ref.dam ref.prof ref.reads.prof ref.minimax.kpos
  ../himapper -v -T1 ref.dam ref.minimax.kpos reads.dam all
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
