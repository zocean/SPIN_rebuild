#!/home/yangz6/Software/Python-2.7.5/python-2.7.5
# Programmer : Yang Zhang 
# Contact: yzhan116@illinois.edu
# Last-modified: 29 Mar 2023 07:49:10 PM

import os,sys,argparse
import pandas as pd

def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--bin',type=str,dest="bin",help="input bin file")
    p.add_argument('--data',type=str,dest="data",help="input data file")
    p.add_argument('--cell_label',type=str,dest="cell_label",nargs="+",help="cell labels match to the input data file")
    p.add_argument('--state',type=str,dest="state",help="state output file")
    p.add_argument('--genome_size',type=str,dest="genome_size",help="genome size file")
    p.add_argument('--output',type=str,dest="output",help="output folder name")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()

def merge_bin(filename):
    with open(filename + '.tmp', 'w') as fout:
        tmp_chrom = None
        tmp_end = None
        tmp_state = None
        tmp_col = None
        with open(filename, 'r') as fin:
            for line in fin:
                if line.strip() == '':
                    continue
                row = line.strip().split()
                chrom = row[0]
                start = row[1]
                end = row[2]
                state = row[3]
                other_col = '\t'.join(row[4:])
                if tmp_chrom is None:
                    tmp_chrom = chrom
                    tmp_start = start
                    tmp_end = end
                    tmp_state = state
                    tmp_col = other_col
                    continue
                if (tmp_chrom != chrom) or (tmp_end != start) or (tmp_state != state):
                    print("%s\t%s\t%s\t%s\t%s" % (tmp_chrom, tmp_start, tmp_end, tmp_state, tmp_col), file = fout)
                    tmp_chrom = chrom
                    tmp_start = start
                    tmp_end = end
                    tmp_state = state
                    tmp_col = other_col
                elif tmp_end == start:
                    tmp_end = end
                else:
                    print("Uknown condition", file = sys.stderr)
        print("%s\t%s\t%s\t%s\t%s" % (tmp_chrom, tmp_start, tmp_end, tmp_state, tmp_col), file = fout)
    os.system("mv %s.tmp %s" % (filename, filename))
 
def bed_to_bigbed(filename, genome_size, bed_type):
    os.system("sort -k1,1 -k2,2n %s >%s.sort" % (filename, filename))
    # merge adjacent bin
    merge_bin("%s.sort"  % (filename))
    os.system("bedToBigBed -type=%s %s.sort %s %s" %(bed_type, filename, genome_size, filename.replace('bed', 'bb')))
    os.unlink("%s.sort" % (filename))
    #os.unlink("%s" % (filename))

def bedgraph_to_bigwig(filename, genome_size):
    os.system("sort -k1,1 -k2,2n %s >%s.sort" % (filename, filename))
    os.system("bedGraphToBigWig %s.sort %s %s" % (filename, genome_size, filename.replace('bedgraph', 'bw')))
    os.unlink("%s.sort" % (filename))
    #os.unlink("%s" % (filename))

def main():
    global args
    args = parse_arg()
    # check output folder 
    if os.path.exists(args.output):
        print("output folder %s already exsits, overwrite content inside" % (args.output), file = sys.stdout)
    else:
        os.mkdir(args.output)
        pass
    # parse result
    bin_list = pd.read_csv(args.bin, sep = '\t', header = None, names = ['chrom', 'start', 'end', 'bin_idx', 'cell'])
    bin_list = bin_list.loc[bin_list['cell'].isin(args.cell_label)]
    state_list = pd.read_csv(args.state, sep = '\t', header = 0, names = ['chrom', 'start', 'end', 'state', 'bin_idx', 'cell'])
    data_list = pd.read_csv(args.data, sep = '\t', header = 0)
    print("Parse data done", file = sys.stdout)
    print("", file = sys.stdout)
    assert bin_list.shape[0] == data_list.shape[0]
    bin_data = pd.concat([bin_list, data_list], axis = 1)
    # merge state with bin
    bin_state = pd.merge(bin_data, state_list, on = ['chrom', 'start', 'end', 'bin_idx', 'cell'], how = 'inner')
    # sort the state by TSA_SON values
    state_TSA_mean = bin_state.groupby(['cell', 'state'])['TSA_SON'].mean().reset_index()
    state_order = list(state_TSA_mean.sort_values(by = 'TSA_SON', ascending = False)['state'])
    print("Reorder state by TSA-seq from high to low")
    print(state_order)
    color_list = ['158,1,66', '220,73,76', '248,141,81', '253,211,128', '255,255,191', '215,239,155', '136,207,164', '63,150,183', '94,79,162']
    state2color = {}
    for state, color in zip(state_order, color_list):
        state2color[str(state)] = color
    # report state
    cell_list = list(state_list['cell'].unique())
    for cell in cell_list:
        data = state_list.loc[state_list['cell'] == cell]
        out_filename = args.output + '/%s_state.bed' % (cell)
        with open(out_filename, 'w') as fout:
            for idx, row in data.iterrows():
                state_label = str(row['state']) 
                print("%s\t%s\t%s\tState_%s\t1000\t+\t%s\t%s\t%s" % (row['chrom'], str(row['start']), str(row['end']), state_label, str(row['start']), str(row['end']), state2color[state_label]), file = fout)
            # sort and convert to bigbed
        bed_to_bigbed(out_filename, args.genome_size, 'bed9')
    # report input data
    cell_list = list(bin_data['cell'].unique())
    label_list = list(data_list.columns)
    for cell in cell_list:
        print("Report input data of cell %s" % (cell), file = sys.stdout)
        data = bin_data.loc[bin_data['cell'] == cell]
        for label in label_list:
            out_filename = args.output + '/%s_%s.bedgraph' % (cell, label)
            data_label = data[['chrom', 'start', 'end', label]]
            data_label.to_csv(out_filename, sep = '\t', header = False, index = False)
            # sort and convert to bigwig
            bedgraph_to_bigwig(out_filename, args.genome_size)
            print("%s done" %(label), file = sys.stdout)
        print("", file = sys.stdout)
 
if __name__=="__main__":
    main()
