#!/bin/bash

# dump Hi-C file
bash scripts/dump -h data/K562_hg38_4DNFITUOMFUQ.hic -j /home/yuchuanw/0_Software/HiC/juicer_tools.1.8.9_jcuda.0.8.jar -l K562 -r 25000 -c hg38_chrom_list -o result/Hi-C

# convert pair to matrix file
bash scripts/pair2matrix -s scripts/pair2matrix.py -b bin_25kb.noblacklist.genome_bin.cell -c K562 -o result/Hi-C/K562/res_25000

# compute p-value 
bash scripts/matrix2pval -s scripts -o result/Hi-C/K562/res_25000

# build final edge for SPIN
python scripts/get_edge.py --chr hg38_chrom_list --bin bin_25kb.noblacklist.genome_bin.cell --cell K562 --folder result/Hi-C/K562/res_25000 --output result/K562.edge

# plot Hi-C input
Rscript scripts/plot_HiC_input.R --bin bin_25kb.noblacklist.genome_bin.cell --cell K562 --hic demo_data/K562_chr1_chr1_oe_VC_SQRT_25000.pair --edge demo_data/K562.edge --chrom_1 chr1 --start_1 40000000 --end_1 50000000 --chrom_2 chr1 --start_2 40000000 --end_2 50000000 --hic_merge_factor 1 --pdf_height 7 --pdf_width 15 demo_data/demo_plot.pdf

# run SPIN

