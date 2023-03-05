#!/home/yangz6/Software/Python-2.7.5/python-2.7.5
# Programmer : Yang Zhang 
# Contact: yzhan116@illinois.edu
# Last-modified: 22 Feb 2023 01:15:17 AM

import os,sys,argparse
import numpy as np

def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--hic',type=str,dest="hic",help="Hi-C O/E or observed result dumped by juicer tool")
    p.add_argument('--bin',type=str,dest="bin",help="bin annotatin file")
    p.add_argument('--cell',type=str,dest="cell",help="cell label")
    p.add_argument('--res',type=int,dest="res",help="resolution")
    p.add_argument('--output',type=str,dest="output",help="output file")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()

def parse_bin_annoation(filename, cell):
    chrom_pos = {}
    pos2id = {}
    with open(filename, 'r') as fin:
        for line in fin:
            row = line.strip().split()
            chrom = row[0]
            start = row[1]
            end = row[2]
            bin_id = int(row[3])
            if row[4] == cell:
                if chrom_pos.get(chrom, None) is None:
                    chrom_pos[chrom] = {'start': 1e10, 'end': 0}
                    pos2id[chrom] = {}
                pos2id[chrom][start] = bin_id
                #
                if chrom_pos[chrom]['start'] > bin_id:
                    chrom_pos[chrom]['start'] = bin_id
                if chrom_pos[chrom]['end'] < bin_id:
                    chrom_pos[chrom]['end'] = bin_id
    return chrom_pos, pos2id

def main():
    global args
    args = parse_arg()
    # parse chrom
    bin_anno,pos2id = parse_bin_annoation(args.bin, args.cell)
    # get chromosome names
    # K562_chr22_chr22_oe_VC_SQRT_25000.matrix
    chr_1 = None
    chr_2 = None
    for item in args.hic.split('/')[-1].split('_'):
        if 'chr' in item:
            if chr_1 is None:
                chr_1 = item
            elif chr_2 is None:
                chr_2 = item
            else:
                pass
    # 
    chr_1_start = bin_anno[chr_1]['start']
    chr_1_end = bin_anno[chr_1]['end']
    chr_2_start = bin_anno[chr_2]['start']
    chr_2_end = bin_anno[chr_2]['end'] 
    result = {}
    with open(args.hic, 'r') as fin:
        for line in fin:
            data = line.strip().split('\t')
            pos_1 = pos2id[chr_1].get(data[0], None)
            pos_2 = pos2id[chr_2].get(data[1], None)
            if pos_1 is None or pos_2 is None:
                continue
            score = float(data[2])
            if np.isnan(score):
                continue
            result[(pos_1, pos_2)] = score
            if chr_1 == chr_2: # intra
                result[(pos_2, pos_1)] = score
    # build
    with open(args.output, 'w') as fout:
        for ii in range(chr_1_start, chr_1_end + 1):
            array = []
            for jj in range(chr_2_start, chr_2_end + 1):
                score = result.get((ii,jj), None)
                if score is None:
                    array.append('0')
                else:
                    array.append('{:.6g}'.format(score))
            print('\t'.join(array), file = fout)

if __name__=="__main__":
    main()
