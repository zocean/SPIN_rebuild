#!/home/yangz6/Software/Python-2.7.5/python-2.7.5
# Programmer : Yang Zhang 
# Contact: yzhan116@illinois.edu
# Last-modified: 26 Feb 2023 09:40:15 PM

import os,sys,argparse
import math
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool

def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--chr',type=str,dest="chr",help="chr list file")
    p.add_argument('--bin_anno',type=str,dest="bin_anno",help="bin annotation file")
    p.add_argument('--cell',type=str,dest="cell",help="cell name")
    p.add_argument('--folder',type=str,dest="folder",help="folder containing the p-value result")
    p.add_argument('--pval_file_suffix',type=str,dest="pval_suffix",default="VC_SQRT_25000.matrix.pvalue",help="suffix of the pvalue file, eg, _VC_SQRT_25000.matrix.pvalue")
    p.add_argument('--output',type=str,dest="output",help="output file")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()

def parse_chr_list(filename):
    # load chromosome list
    chr_list = []
    with open(filename, 'r') as fin:
        for line in fin:
            if line.strip() == '':
                continue
            else:
                chr_list.append(line.strip().split()[0])
    return chr_list

def get_bin_idx_per_chrom(filename, cell_id):
    if cell_id is None:
        print("Cell name must be provided", file = sys.stderr)
        exit(1)
    # get the range of bin_idx for each chromosome
    bin_idx_table = {}
    idx2pos = {}
    with open(filename, 'r') as fin:
        for line in fin:
            if line.strip() == '':
                continue
            row = line.strip().split()
            chrom = row[0]
            start= row[1]
            end = row[2]
            idx = int(row[3])
            cell = row[4]
            idx2pos[idx] = (chrom, start, end)
            if cell == args.cell:
                #
                if bin_idx_table.get(chrom, None) is None:
                    bin_idx_table[chrom] = [idx,-1] # set the first index
                else:
                    #
                    if idx < bin_idx_table[chrom][0]: # 
                        print("idx range error", file = sys.stderr)
                    bin_idx_table[chrom][1] = idx # update the index to get the larget one
            else:
                pass
    return bin_idx_table, idx2pos

def smooth_matrix(matrix, pad_width = 1, method = 'mean'):
    # pad matriix
    padded_matrix = np.pad(matrix, pad_width, 'constant', constant_values = np.nan)
    window_size = pad_width + 2
    window = np.lib.stride_tricks.sliding_window_view(padded_matrix, (window_size, window_size)).reshape((matrix.shape[0], matrix.shape[1], -1))

    if method == 'median':
        # compute the median
        result = np.nanmedian(window, axis = 2)
    elif method == 'mean':
        # compute the mean
        result = np.nanmean(window, axis = 2)
    else:
        result = np.nanmean(window, axis = 2)

    return result
    
def extract_edge(pvalue_file, chr_idx_table, idx2pos, chr_i, chr_j, is_smooth, smooth_method, pvalue_cutoff):
    edge_list = []
    # load data as matrix
    pval_matrix = np.loadtxt(pvalue_file)
    # check the dimension of the matrix
    true_dim_i = chr_idx_table[chr_i][1] - chr_idx_table[chr_i][0] + 1
    true_dim_j = chr_idx_table[chr_j][1] - chr_idx_table[chr_j][0] + 1
    try:
        assert true_dim_i == pval_matrix.shape[0] and true_dim_j == pval_matrix.shape[1]
    except AssertionError:
        print("Pvalue matrix result dimension does not match", file = sys.stderr) 
        print(pval_matrix.shape)
        print(true_dim_i, true_dim_j)
    # smooth if needed
    if is_smooth:
        final_matrix = smooth_matrix(pval_matrix, pad_width = 1, method = smooth_method) # smooth by first neighbor
    else:
        final_matrix = pval_matrix
    try:
        assert final_matrix.shape == pval_matrix.shape
    except AssertionError:
        print("Matrix dimension is wrong", file = sys.stderr)
    # extract final edges
    for ii in range(true_dim_i):
        for jj in range(true_dim_j):
            bin_idx_i = chr_idx_table[chr_i][0] + ii
            bin_idx_j = chr_idx_table[chr_j][0] + jj
            log_pval = final_matrix[ii, jj]
            # do nothing on diaganol
            if chr_i == chr_j and ii == jj:
                continue
            #  connect adjacent bins unless the original annotation says they are not continous on the genome
            elif chr_i == chr_j and abs(ii-jj) < 1.9:
                if ii < jj and idx2pos[ii][2] == idx2pos[jj][1]:
                    edge_list.append((bin_idx_i, bin_idx_j, pvalue_cutoff))
                elif jj < ii and idx2pos[jj][2] == idx2pos[ii][1]:
                    edge_list.append((bin_idx_i, bin_idx_j, pvalue_cutoff))
                else:
                    pass
            #
            else:
                try:
                    if log_pval >= pvalue_cutoff:
                        edge_list.append((bin_idx_i, bin_idx_j, log_pval))
                except ValueError:
                    print(log_pval)
    return chr_i, chr_j, edge_list

def main():
    global args
    args = parse_arg()
    # 
    chr_list = parse_chr_list(args.chr) 
    # 
    chr_idx_table,idx2pos = get_bin_idx_per_chrom(args.bin_anno, args.cell)
    # main function to filter
    N = len(chr_list)
    pvalue_cutoff = -1 * math.log(1e-5)
    print("-log(pvalue) cutoff is %s" % (pvalue_cutoff), file = sys.stdout)

    progress_bar = tqdm(total= (N*(N-1)/2 + N))

    result_table = {}
     
    def log_result(result):
        result_table[(result[0], result[1])] = result[2]
        progress_bar.update(1)

    def error_result(error):
        print("Get error: %s" % (error), flush = True)

    pool = Pool(processes = 16)
    for ii in range(N):
        for jj in range(N)[ii:]:
            chr_i = chr_list[ii]
            chr_j = chr_list[jj]
            if chr_i == chr_j:
                pvalue_file = os.path.join(args.folder, "%s_%s_%s_oe_%s" % (args.cell, chr_i, chr_j, args.pval_suffix))
            else:
                pvalue_file = os.path.join(args.folder, "%s_%s_%s_observed_%s" % (args.cell, chr_i, chr_j, args.pval_suffix))
            if not os.path.isfile(pvalue_file):
                print("Can't open file: %s" % (pvalue_file), file = sys.stderr)
                continue
            # intra-chromosomal interactions
            if chr_i == chr_j: 
                is_smooth = False
            else: # inter-chromosomomal interactions
                is_smooth = True
            pool.apply_async(extract_edge, [pvalue_file, chr_idx_table, idx2pos, chr_i, chr_j, is_smooth, 'mean', pvalue_cutoff], callback = log_result, error_callback = error_result)
    pool.close()
    pool.join()
 
    # report to output file
    with open(args.output, 'w') as fout:
        for key in sorted(result_table.keys()):
            edge_list = result_table[key]
            for edge in edge_list:
                print("%d\t%d\t%.4f" % (edge[0], edge[1], edge[2]), file = fout)
        
if __name__=="__main__":
    main()
