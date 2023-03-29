#!/bin/bash


# merge all input
head -n 1 result/K562.input.imputed >result/4cell.input.imputed.header
cat result/K562.input.imputed result/HCT116.input.imputed result/H1.input.imputed result/HFF.input.imputed | grep -v TSA >result/4cell.input.imputed
cat result/4cell.input.imputed.header result/4cell.input.imputed >result/4cell.input.imputed.merge
rm result/4cell.input.imputed.header
mv result/4cell.input.imputed.merge result/4cell.input.imputed

# baseline hmm
python scripts/baseline_hmm.py --input result/4cell.input.imputed --bin bin_25kb.noblacklist.genome_bin.cell --n_state 9 --seed 10 --merge_cutoff 250000 --cell_label K562 HCT116 H1 HFF  --output result/4cell.state.hmm


