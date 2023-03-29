#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 29 Mar 2023 06:33:33 PM

import os,sys,argparse
import math

def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--script',type=str,dest="script",help="plotting script")
    p.add_argument('--bin_anno',type=str,dest="bin_anno",help="bin annotation file")
    p.add_argument('--cell_label',type=str,dest="cell_label",help="cell label to used")
    p.add_argument('--edge',type=str,dest="edge",help="edge file")
    p.add_argument('--hic',type=str,dest="hic",help="folder contains the hic pair file")
    p.add_argument('--hic_merge_factor',type=int,dest="hic_merge_factor",help="pooling factor")
    p.add_argument('--genome',type=str,dest="genome",help="chromosome size file")
    p.add_argument('--output',type=str,dest="output",help="output folder")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()

def parse_genome(filename):
    genome_table = {}
    with open(filename, 'r') as fin:
        for line in fin:
            if line.strip() == '':
                continue
            row = line.strip().split()
            chrom = row[0]
            size = int(row[1])
            genome_table[chrom] = size
    return genome_table

def chrom_filter(chrom):
    """
    remove chromosome not in autosomes and sex chromosomes
    """
    if 'Un' in chrom or 'random' in chrom or 'hap' in chrom or 'chrM' in chrom:
        return False
    else:
        return True

def main():
    global args
    args = parse_arg()
    # parse the chromosome size file
    genome_table = parse_genome(args.genome)
    # get chromosome list
    chrom_list= [chrom for chrom in sorted(genome_table.keys()) if chrom_filter(chrom)]
    print("# Generate intra-chromosomal contact map for chromosome:", file = sys.stdout)
    print("# " + ', '. join(chrom_list), file = sys.stdout)
    # build script
    for chrom in chrom_list:
        # K562_chr2_chr21_observed_VC_SQRT_25000.pair
        pair_file = "%s/%s_%s_%s_oe_VC_SQRT_25000.pair" % (args.hic, args.cell_label, chrom, chrom)
        if not os.path.isfile(pair_file):
            print("Can't find pair file %s" % (pair_file), file = sys.stderr)
            continue
        #
        start = 0
        end = int(math.floor(genome_table[chrom] / 1000000)*1000000)
        #
        pdf_width = int(round( (end - start) / 1000000 / 5) + 1)
        pdf_height = int(round( (end - start) / 1000000 / 10))
        cmd = "Rscript %s --bin %s --cell %s --hic %s --edge %s --chrom_1 %s --start_1 %d --end_1 %d --chrom_2 %s --start_2 %d --end_2 %d --hic_merge_factor %d --pdf_height %d --pdf_width %d %s" % (args.script, args.bin_anno, args.cell_label, pair_file, args.edge, chrom, start, end, chrom, start, end, args.hic_merge_factor, pdf_height, pdf_width, os.path.join(args.output, "%s_%s-%s.pdf" % (args.cell_label, chrom, chrom)))
        print(cmd, file = sys.stdout)

if __name__=="__main__":
    main()
