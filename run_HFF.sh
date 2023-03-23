#!/bin/bash

# input parameters
HIC_FILE=data/H1_LMNB1_DAMID_25kb_hg38_combined.bw
CELL=HFF
TSA_SON=U54_data/hff/TSA_Seq/25kb/HFFc6_SON_TSA_25kb_hg38_rep1_201806condE.bw
TSA_LMNB1=U54_data/hff/TSA_Seq/25kb/HFFc6_LMNB1_TSA_25kb_hg38_rep1_201810condAI.bw
TSA_MKI67IP=U54_data/hff/TSA_Seq/25kb/HFFc6_LMNB1_TSA_25kb_hg38_rep1_201810condAI.bw
DAMID_LMNB1=U54_data/hff/DamID/25kb/HFFc6_LMNB1_TSA_25kb_hg38_rep1_201810condAI.bw

# dump Hi-C file
bash scripts/dump -h $HIC_FILE -j /home/yuchuanw/0_Software/HiC/juicer_tools.1.8.9_jcuda.0.8.jar -l $CELL -r 25000 -c hg38_chrom_list -o result/Hi-C

# convert pair to matrix file
bash scripts/pair2matrix -s scripts/pair2matrix.py -b bin_25kb.noblacklist.genome_bin.cell -c $CELL -o result/Hi-C/$CELL/res_25000

# compute p-value 
bash scripts/matrix2pval -s scripts -o result/Hi-C/$CELL/res_25000

# build final edge for SPIN
python scripts/get_edge.py --chr hg38_chrom_list --bin bin_25kb.noblacklist.genome_bin.cell --cell $CELL --folder result/Hi-C/$CELL/res_25000 --output result/$CELL.edge

# plot Hi-C input
#Rscript scripts/plot_HiC_input.R --bin bin_25kb.noblacklist.genome_bin.cell --cell K562 --hic demo_data/K562_chr1_chr1_oe_VC_SQRT_25000.pair --edge demo_data/K562.edge --chrom_1 chr1 --start_1 40000000 --end_1 50000000 --chrom_2 chr1 --start_2 40000000 --end_2 50000000 --hic_merge_factor 1 --pdf_height 7 --pdf_width 15 demo_data/demo_plot.pdf

# process the TSA-seq data
python scripts/power_transform.py --bw $TSA_Seq $TSA_LMNB1 $TSA_MKI67IP $DAMID_LMNB1 --label TSA_SON TSA_LMNB1 TSA_MKI67IP DamID_LMNB1 --bin bin_25kb.noblacklist.genome_bin.cell --cell_name $CELL --output result/$CELL.input --original result/$CELL.input.original

# impurte missing input data
Rscript scripts/impute_NA.R --output result/$CELL.input.imputed --na_string NA --seed 11 result/$CELL.input 2>&1 >result/$CELL.input.imputed.log

# run HMM baseline 
python scripts/baseline_hmm.py --input result/$CELL.input.imputed --bin bin_25kb.noblacklist.genome_bin.cell --n_state 9 --seed 10 --merge_cutoff 250000 --cell_label $CELL --output result/$CELL.state.hmm 
python scripts/baseline_hmm.py --input result/$CELL.input.imputed --bin bin_25kb.noblacklist.genome_bin.cell --n_state 9 --seed 11 --merge_cutoff 250000 --cell_label $CELL --output result/$CELL.state.hmm.seed11
python scripts/baseline_hmm.py --input result/$CELL.input.imputed --bin bin_25kb.noblacklist.genome_bin.cell --n_state 9 --seed 12 --merge_cutoff 250000 --cell_label $CELL --output result/$CELL.state.hmm.seed12

# plot HMM result 
#$Rscript scripts/plot_HiC_state.R --bin bin_25kb.noblacklist.genome_bin.cell --cell $CELL --input result/K562.input.imputed --state result/K562.state.hmm --edge demo_data/K562.edge --chrom_1 chr1 --start_1 40000000 --end_1 50000000 --chrom_2 chr1 --start_2 40000000 --end_2 50000000 --hic_merge_factor 1 --pdf_height 7 --pdf_width 15 demo_data/demo_plot_hmm.pdf

# run SPIN
