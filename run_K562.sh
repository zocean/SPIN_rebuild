#!/bin/bash

# input parameters
HIC_FILE=data/K562_hg38_4DNFITUOMFUQ.hic
CELL=K562
TSA_SON=U54_data/k562/TSA_Seq/25kb/K562_SON_TSA_25kb_hg38_rep1_201911condE.bw
TSA_LMNB1=U54_data/k562/TSA_Seq/25kb/K562_LMNB_TSA_25kb_hg38_rep1_201712condA.bw
TSA_MKI67IP=U54_data/k562/TSA_Seq/25kb/K562_WT_MKI67IP_2020_08_rep1_25kb_hg38.bw
DAMID_LMNB1=U54_data/k562/DamID/25kb/K562_LMNB1_DAMID_25kb_hg38_combined.bw

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
Rscript scripts/plot_HiC_input.R --bin bin_25kb.noblacklist.genome_bin.cell --cell K562 --hic demo_data/K562_chr2_chr2_oe_VC_SQRT_25000.pair --edge demo_data/K562.edge --chrom_1 chr2 --start_1 20000000 --end_1 90000000 --chrom_2 chr2 --start_2 20000000 --end_2 90000000 --hic_merge_factor 8 --pdf_height 7 --pdf_width 15 demo_data/demo_plot2.pdf

# process the TSA-seq data
python scripts/power_transform.py --bw $TSA_Seq $TSA_LMNB1 $TSA_MKI67IP $DAMID_LMNB1 --label TSA_SON TSA_LMNB1 TSA_MKI67IP DamID_LMNB1 --bin bin_25kb.noblacklist.genome_bin.cell --cell_name $CELL --output result/$CELL.input --original result/$CELL.input.original

# impurte missing input data
Rscript scripts/impute_NA.R --output result/$CELL.input.imputed --na_string NA --seed 11 result/$CELL.input 2>&1 >result/$CELL.input.imputed.log

# run HMM baseline 
python scripts/baseline_hmm.py --input result/$CELL.input.imputed --bin bin_25kb.noblacklist.genome_bin.cell --n_state 9 --seed 10 --merge_cutoff 250000 --cell_label $CELL --output result/$CELL.state.hmm 
python scripts/baseline_hmm.py --input result/$CELL.input.imputed --bin bin_25kb.noblacklist.genome_bin.cell --n_state 9 --seed 11 --merge_cutoff 250000 --cell_label $CELL --output result/$CELL.state.hmm.seed11
python scripts/baseline_hmm.py --input result/$CELL.input.imputed --bin bin_25kb.noblacklist.genome_bin.cell --n_state 9 --seed 12 --merge_cutoff 250000 --cell_label $CELL --output result/$CELL.state.hmm.seed12

# create trasck
python scripts/build_session.py --bin bin_25kb.noblacklist.genome_bin.cell --cell_label K562 --data result/K562.input.imputed --state result/K562.state.hmm --genome hg38.genome --output result/browser_track/ 

# plot HMM result 
#$Rscript scripts/plot_HiC_state.R --bin bin_25kb.noblacklist.genome_bin.cell --cell $CELL --input result/K562.input.imputed --state result/K562.state.hmm --edge demo_data/K562.edge --chrom_1 chr1 --start_1 40000000 --end_1 50000000 --chrom_2 chr1 --start_2 40000000 --end_2 50000000 --hic_merge_factor 1 --pdf_height 7 --pdf_width 15 demo_data/demo_plot_hmm.pdf

# run SPIN

# run annoation
# RT
grep 'K562' bin_25kb.noblacklist.genome_bin.cell | python scripts/bed_anno.py --bed - --anno ../data/RT_multi-fraction/K562_7frac_repli-seq_rep1_hg38_G2.bw ../data/RT_multi-fraction/K562_7frac_repli-seq_rep1_hg38_S1.bw ../data/RT_multi-fraction/K562_7frac_repli-seq_rep1_hg38_S2.bw ../data/RT_multi-fraction/K562_7frac_repli-seq_rep1_hg38_S3.bw ../data/RT_multi-fraction/K562_7frac_repli-seq_rep1_hg38_S4.bw ../data/RT_multi-fraction/K562_7frac_repli-seq_rep1_hg38_S5.bw ../data/RT_multi-fraction/K562_7frac_repli-seq_rep1_hg38_S6.bw --mode signal_mean signal_mean signal_mean signal_mean signal_mean signal_mean signal_mean --label G2 S1 S2 S3 S4 S5 S6 --genome ../annotation/hg38.fa --genome_size ../annotation/hg38.genome --output result/annotation/K562_multi-frac_RT.tsv
# histone

grep 'K562' bin_25kb.noblacklist.genome_bin.cell | python scripts/bed_anno.py --bed - --anno ../data/Histone/K562_ChIP-seq_H2AFZ_foldChange_GRCh38_ENCFF621DJP.bw ../data/Histone/K562_ChIP-seq_H3K27ac_foldChange_GRCh38_ENCFF849TDM.bw ../data/Histone/K562_ChIP-seq_H3K27me3_foldChange_GRCh38_ENCFF405HIO.bw ../data/Histone/K562_ChIP-seq_H3K36me3_foldChange_GRCh38_ENCFF163NTH.bw ../data/Histone/K562_ChIP-seq_H3K4me1_foldChange_GRCh38_ENCFF834SEY.bw ../data/Histone/K562_ChIP-seq_H3K4me2_foldChange_GRCh38_ENCFF959YJV.bw ../data/Histone/K562_ChIP-seq_H3K4me3_foldChange_GRCh38_ENCFF660WUG.bw ../data/Histone/K562_ChIP-seq_H3K79me2_foldChange_GRCh38_ENCFF957YJT.bw ../data/Histone/K562_ChIP-seq_H3K9ac_foldChange_GRCh38_ENCFF286WRJ.bw ../data/Histone/K562_ChIP-seq_H3K9me3_foldChange_GRCh38_ENCFF601JGK.bw ../data/Histone/K562_ChIP-seq_H4K20me1_foldChange_GRCh38_ENCFF605FAF.bw --mode signal_mean signal_mean signal_mean signal_mean signal_mean signal_mean signal_mean signal_mean signal_mean signal_mean signal_mean --label H2AZ H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me2 H3K4me3 H3K79me2 H3K9ac H3K9me3 H4K20me1 --genome ../annotation/hg38.fa --genome_size ../annotation/hg38.genome --output result/annotation/K562_histone.tsv



