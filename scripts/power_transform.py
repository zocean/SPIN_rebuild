#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 09 Mar 2023 10:03:08 PM

import os,sys,argparse
from bx.bbi.bigwig_file import BigWigFile
from tqdm import tqdm
import numpy as np
from sklearn.preprocessing import PowerTransformer

def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--bw',type=str,dest="bw",nargs="+",help="bigwig file list")
    p.add_argument('--label',type=str,dest="label",nargs="+",help="label for the bigwig files")
    p.add_argument('--bin',type=str,dest="bin",help="bin annotation file")
    p.add_argument('--cell_name',type=str,dest="cell_name",help="cell name to extract in the bin annotation file, only rows in the fifth column of the bin annotation matching the cell_name will be extracted")
    p.add_argument('--original',type=str,dest="original",help="original signal without normalization")
    p.add_argument('--output',type=str,dest="output",help="output file")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()

class Bin(object):
    def __init__(self, chrom, start, end, idx):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.bin_idx = idx
        #
        self.anno = {}
        self.normalized_anno = {}
    def __repr__(self):
        return "%s\t%d\t%d\t%s" % (self.chrom, self.start, self.end, self.bin_idx)
    def cal_mean_signal(self, bw, label):
        '''get the mean signal value within region'''
        array = np.array([])
        signal = bw.get_as_array(bytes(self.chrom, 'utf-8'), self.start, self.end)
        if signal is not None:
            array = np.append(array, signal)
        if len(array) < 1: # resolve the warning issue with numpy when array with zero element
            value = np.NAN
        else:
            value = np.nanmean(array)
        #
        self.anno[label] = value

def parse_bin_anno(filename, cell_name):
    region_list = []
    with open(filename, 'r') as fin:
        for line in fin:
            if line.strip().startswith('#') or line.strip() == '':
                continue
            row = line.strip().split('\t')
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            bin_idx = row[3] 
            cell_label = row[4]
            if cell_label == cell_name:
                region_list.append(Bin(chrom, start, end, bin_idx))
    return region_list

def main():
    global args
    args = parse_arg()
    # check parameters
    try:
        assert len(args.bw) == len(args.label)
    except AssertionError:
        print("length of bigwig files and label does not match", file = sys.stderr)
        exit(1)
    bw_list = []
    label_list = []
    for nn in range(len(args.bw)):
        bw_list.append(BigWigFile(open(args.bw[nn], 'rb')))
        label_list.append(args.label[nn])
    # parse bin annotation 
    region_list = parse_bin_anno(args.bin, args.cell_name)
    print("Parse %d region" % (len(region_list)), file = sys.stdout)
    # annotate bin with bigwig files
    for nn in tqdm(range(len(region_list))):
        region = region_list[nn]
        for bw,label in zip(bw_list, label_list):
            region.cal_mean_signal(bw, label)
    print("calculate mean signal per bin done", file = sys.stdout)
    # power transformation
    for label in label_list:
        value_list = []
        idx_list = []
        for idx in tqdm(range(len(region_list))):
            idx_list.append(idx)
            region = region_list[idx]
            value = region.anno.get(label, np.NAN)
            value_list.append(value)
        value_list = np.array(value_list).reshape(-1, 1)
        pt = PowerTransformer(method='yeo-johnson', standardize=True)
        value_list_transformed = pt.fit_transform(value_list)
        assert value_list.shape[0] == value_list_transformed.shape[0]
        for idx,transformed_value in zip(idx_list, value_list_transformed.flatten()):
            region_list[idx].normalized_anno[label] = transformed_value
    print("transformation done", file = sys.stdout)
    # report 
    if args.original is not None:
        with open(args.original, 'w') as fout:
            print('\t'.join(label_list), file = fout)
            for region in region_list:
                array = []
                for label  in label_list:
                    value = region.anno[label]
                    if np.isnan(value):
                        value = 'NA'
                    else:
                        value = '{:.6f}'.format(value)
                    array.append(value)
                print('\t'.join(array), file = fout)
    with open(args.output, 'w') as fout:
        print('\t'.join(label_list), file = fout)
        for region in region_list: 
            array = []
            for label in label_list:
                value = region.normalized_anno[label]
                if np.isnan(value):
                    value = 'NA'
                else:
                    value = '{:.6f}'.format(value)
                array.append(value)
            print('\t'.join(array), file = fout)
 
if __name__=="__main__":
    main()
