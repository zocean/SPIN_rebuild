#!/home/yangz6/Software/miniconda3/envs/base2/bin/python
# Programmer : Yang Zhang 
# Contact: yangz6@andrew.cmu.edu
# Last-modified: 09 Mar 2023 10:38:05 PM

import os,sys,argparse
from tqdm import tqdm
import math

def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--bw',type=str,dest="bw",help="bigwig file")
    p.add_argument('--genome',type=str,dest="genome",help="chromosome size file")
    p.add_argument('--smooth_res',type=str,dest="smooth_res",help="smooth resolution in bp")
    p.add_argument('--smooth_win',type=int,dest="smooth_win",help="window size for smoothing the data")
    p.add_argument('--output',type=str,dest="output",help="output file prefix")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()

def parse_genome(genome_file):
    table = {}
    with open(genome_file, 'r') as fin:
        for line in fin:
            row = line.strip().split()
            table[row[0]] = int(row[1])
    return table

def smooth(x,window_len,window='hanning'):
    if x.ndim != 1:
        print("smooth only accepts 1 dimension arrays.", file = sys.stderr)
    if x.size < window_len:
        print("Input vector needs to be bigger than window size.", file = sys.stderr)
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        print("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'", file = sys.stderr)
    s = np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]

def main():
    global args
    args = parse_arg()
    # load genome
    genome_table = parse_genome(genome_file)
    # load original bigwig file
    bw = BigWgiFile(open(bw_file, 'rb'))
    # 
    fout = WriteToFile(args.output + '.bedgraph')
    for chrom in sorted(genome_tables.keys()):
        chrom_size = chr_table[chrom]
        value_list = []
        for nn in tqdm(range(math.floor(chrom_size/args.smooth_res))):
            start = int(nn*smooth_res)
            end = int(start + smooth_res)
            array = bw.get_as_array(bytes(chrom, 'utf-8'), start, end)
            value = np.nanmean(array)
            if np.isnan(value):
                value = 0.0
            value_list.append(value)
        agg_array = np.array(value_list)
        smooth_array = smooth(agg_array, args.smooth_win)
        for nn,value in enumerate(smooth_array):
            if nn == 0: 
                print >>fout, "%s\t0\t%d\t%.6f" % (chrom, (nn+1)*args.smooth_res, float(value))
            elif nn == len(smooth_array) - 1:
                print >>fout, "%s\t%d\t%d\t%.6f" % (chrom,nn*args.smooth_res, chrom_size, float(value))
            else:
                print >>fout, "%s\t%d\t%d\t%.6f" % (chrom,nn*args.smooth_res, (nn+1)*args.smooth_res, float(value))
    fout.close()
    # sort and convert bedgraph to bigwig
    cmd = "sort -k1,1 -k2,2n %s >%s" % (args.output + '.bedgraph', args.output + '.bedgraph.sort')
    os.system(cmd)
    cmd = "wigToBigWig -clip %s %s %" % (args.output + '.bedgraph.sort', args.genome, args.output + '.bw')
    os.system(cmd)
    os.unlink(args.output + '.bedgraph')
    os.unlink(args.output + '.bedgraph.sort')

if __name__=="__main__":
    main()
