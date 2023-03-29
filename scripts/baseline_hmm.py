#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 29 Mar 2023 07:48:25 PM

import os,sys,argparse
import numpy as np
from sklearn.experimental import enable_iterative_imputer
from hmmlearn.hmm import GaussianHMM

def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--input',type=str,dest="input",help="SPIN input file")
    p.add_argument('--bin',type=str,dest="bin",help="bin annotation file, the number of lines of this file should be the same of the input file")
    p.add_argument('--cell_label',type=str,dest="cell_label",nargs="+",help="only bins from these cell(s) will be extracted, the fifth column of --bin contains the cell name")
    p.add_argument('--n_state',type=int,dest="n_state",help="number of state")
    p.add_argument('--seed',type=int,dest="seed",default=11,help="random state seed")
    p.add_argument('--merge_cutoff',type=int,dest="merge_cutoff",help="distance cutoff to merge adjacent bins when learn HMM modal (ie, assuming these bins are continous)")
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
        self.cell = None

def parse_bin(filename, cell_list):
    bin_list = []
    with open(filename, 'r') as fin:
        for line in fin:
            if line.strip().startswith('#') or line.strip() == '':
                continue
            row = line.strip().split()
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])
            bin_idx = row[3]
            cell_label = row[4]
            if cell_label in cell_list:
                region = Bin(chrom, start, end, bin_idx)
                region.cell = cell_label
                bin_list.append(region)
    return bin_list

def parse_input(filename, header = True):
    header_table = {}
    N = 0
    data_list = []
    N_NA = 0
    with open(filename, 'r') as fin:
        for line in fin:
            if line.strip().startswith('#') or line.strip() == '':
                continue
            row = line.strip().split()
            if header and N == 0:
                for nn in range(len(row)):
                    header_table[nn] = row[nn]
                N += 1
                continue
            data = []
            for item in row:
                if item == 'NA' or item == 'nan':
                    value = np.NAN
                else:
                    try:
                        value = float(item)
                    except ValueError:
                        print("Can't conver value to float: %s" % (item), file = sys.stderr)
                        exit(1)
                data.append(value)
            assert len(data) == len(header_table)
            if np.NAN in data:
                N_NA += 1
            data_list.append(data)
    print("%d bins parsed, %d bins have missing data " % (len(data_list), N_NA), file = sys.stderr)
    return np.array(data_list)

def split_data_by_gap(bin_clean, merge_cutoff):
    """
    bin annotion file
    merge_cutoff: max distance to merge non-continus bins
    """
    gap_list = []
    tmp_chrom = None
    tmp_end = None
    count = 0
    for region in bin_clean:
        if tmp_chrom is None:
            tmp_chrom = region.chrom
            tmp_end = region.end
            count = 1
        else:
            if tmp_chrom == region.chrom:
                if tmp_end == region.start:
                    tmp_end = region.end
                    count += 1
                else:
                    distance = region.start - tmp_end
                    if distance < merge_cutoff:
                        tmp_end = region.end
                        count += 1
                    else:
                        gap_list.append(count)
                        tmp_chrom = region.chrom
                        tmp_end = region.end
                        count = 1
            else:
                gap_list.append(count)
                tmp_chrom = region.chrom
                tmp_end = region.end
                count = 1
    gap_list.append(count)
    count = 0
    assert int(np.sum(gap_list)) == len(bin_clean)
    return np.array(gap_list)

def main():
    global args
    args = parse_arg()
    print("Randome seed is %d" % (args.seed), file = sys.stdout)
    np.random.seed(args.seed)
    print("Non-continous bin merge cutoff is %d" % (args.merge_cutoff), file = sys.stdout)
    # parse the bin annotation file
    bin_list = parse_bin(args.bin, args.cell_label)
    # parse the input data
    data = parse_input(args.input)
    # check the number of rows of bin and input data are the same
    try:
        assert len(bin_list) == data.shape[0]
    except AssertionError:
        print("Number of bin annotation file (%d) and input data (%d) do not match" % (len(bin_list), data.shape[0]), file = sys.stderr)
        exit(1)
    # filter rows without NA
    clean_idx = np.arange(data.shape[0])[~np.isnan(data).any(axis = 1)]
    data_clean = data[clean_idx]
    bin_clean = [bin_list[idx] for idx in clean_idx]
    print("%d rows are removed because of missing values" % (data.shape[0] - data_clean.shape[0]), file = sys.stdout)
    # split the data into gapped squence
    gap_list = split_data_by_gap(bin_clean, args.merge_cutoff)
    # fit HMM model
    print("Number of states is %d" % (args.n_state))
    init_startprob = np.ones(args.n_state)/args.n_state
    init_transmat = np.ones([args.n_state, args.n_state])/args.n_state
    model = GaussianHMM(n_components = args.n_state, covariance_type = 'full', n_iter = 100, random_state = args.seed)
    model.startprob_ = init_startprob
    model.transmat_ = init_transmat
    model.fit(data_clean, gap_list)
    result = model.predict(data_clean, gap_list) 
    # 
    assert result.shape[0] == len(bin_clean)
    with open(args.output, 'w') as fout:
        print("chrom\tstart\tend\tstate\tbin\tcell", file = fout)
        for idx in range(len(bin_clean)):
            region = bin_clean[idx]
            state = result[idx]
            print("%s\t%d\t%d\t%d\t%s\t%s" % (region.chrom, region.start, region.end, state, region.bin_idx, region.cell), file = fout)
 
if __name__=="__main__":
    main()
