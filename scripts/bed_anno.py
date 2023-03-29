#!/home/yangz6/Software/Python-2.7.5/python-2.7.5
# Programmer : Yang Zhang 
# Contact: yzhan116@illinois.edu
# Last-modified: 29 Mar 2023 07:32:52 PM

import os,sys,argparse
from tqdm import tqdm
import pysam
from bx.bbi.bigwig_file import BigWigFile
import tabix
import pypairix
import random
'''import custom function/class'''
from utility import *
from multiprocessing import Pool, Manager

def parse_arg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -h', epilog='Library dependency :')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('--thread',type=int,dest="thread",default=16,help="number of threads used")
    p.add_argument('--bed',type=str,dest="bed",help="bed file also support stdin use '-' to indicate the bed file is from stdin")
    p.add_argument('--update',dest="update",action="store_true",help="if set, will update the file assumeing header exists")
    p.add_argument('--exclude',type=str,dest="exclude",nargs="+",help="exclude region in tabix format")
    p.add_argument('--anno',type=str,dest="anno",nargs="+",help="annotation file list")
    p.add_argument('--mode',type=str,dest="mode",nargs="+",help="mode list")
    p.add_argument('--label',type=str,dest="label",nargs="+",help="label list")
    p.add_argument('--genome',type=str,dest="genome",help="genome fasta file")
    p.add_argument('--genome_size',type=str,dest="genome_size",help="genome size file")
    p.add_argument('--resort',dest="resort",action="store_true",help="if set this argument, the output bed fill bed sorted by chrom then by the start position")
    p.add_argument('--output',type=str,dest="output",help="output file name")
    if len(sys.argv) < 2:
        print(p.print_help())
        exit(1)
    return p.parse_args()

def load_exclude(file_list):
    exclude_list = []
    if file_list is None:
        return exclude_list
    for filename in file_list:
        exclude_list.append(tabix.open(filename))
    return exclude_list

def options_check():
    text = ""
    text += "# Bed file: %s\n" % (args.bed)
    if args.update:
        text += "# update existing annotation: True\n"
    else:
        text += "# update existing annotation: False\n"
    text += "# genome fasta file: %s\n" % (args.genome)
    text += "# genome chromosome size file: %s\n" % (args.genome_size)
    if args.exclude is not None:
        text += "# exclude region file: " + ','.join(filename for filename in args.exclude) + '\n'
    text += "# output file: %s\n" % (args.output)
    for nn in range(len(args.anno)):
        anno = args.anno[nn]
        label = args.label[nn]
        mode = args.mode[nn]
        text += "# Annotaion label: %s\n" % (label)
        text += "#\t\tAnnotation mode: %s\n" % (mode)
        text += "#\t\tAnnotation file: %s\n" % (anno)
    print(text, file = sys.stderr)

def single_region_run(region):
    for nn in range(len(anno_list)):
        label = label_list[nn]
        mode = mode_list[nn]
        if mode in ['coverage', 'coverage_per', 'count', 'anno', 'nearby']:
            anno = tabix.open(anno_list[nn])
        elif mode in ['loop']:
            anno = pypairix.open(anno_list[nn])
        elif mode in ['signal_mean', 'signal_sum']:
            anno = BigWigFile(open(anno_list[nn], 'rb'))
        elif mode in ['dnase_count','dnase_ave']:
            anno = SamHandler(pysam.Samfile(anno_list[nn]))
        elif mode in ['gc', 'gc_count']:
            anno = pysam.Fastafile(anno_list[nn])
        elif mode in ['size']:
            anno = anno_list[nn]
        else:
            print("Unknown mode: %s" % (mode))
        try:
            region.get_anno(anno, label, mode, genome_size)
        except:
            print("[anno] skip %s:%d-%d" % (region.chrom, region.start, region.end), file = sys.stderr)
            continue
    return region

def Main():
    global args
    args = parse_arg() # check parameter
    if args.anno is None and args.label is None and args.mode is None:
        args.anno = []; args.label = []; args.mode = []
        pass
    else:
        assert len(args.anno) == len(args.label)
        assert len(args.anno) == len(args.mode)
    options_check()
    global anno_list, label_list, mode_list, genome_size
    # load genome fasta/size
    genome = pysam.Fastafile(args.genome)
    genome_size = load_genome_size(args.genome_size)
    # load exclude tabix list
    if args.exclude is not None:
        print("load exclude regions", file = sys.stderr)
        exclude_list = load_exclude(args.exclude)
    else:
        exclude_list = []
    # load bed region
    if args.update:
        region_list = load_region_anno(args.bed)
    else:
        region_list = load_region(args.bed)
    # load annotation files into tabix, bigwig, fasta or bam object
    anno_list = []
    label_list = args.label
    mode_list = args.mode
    anno_list = args.anno
    #for nn in range(len(mode_list)):
    #    mode = mode_list[nn]
    #    if mode in ['coverage', 'coverage_per', 'count', 'anno', 'nearby']:
    #        try:
    #            anno_list.append(tabix.open(args.anno[nn]))
    #        except:
    #            print(args.anno[nn], file = sys.stderr)
    #            exit(1)
    #    elif mode in ['loop']:
    #        anno_list.append(pypairix.open(args.anno[nn]))
    #    elif mode in ['signal_mean', 'signal_sum']:
    #        anno_list.append(BigWigFile(open(args.anno[nn], 'rb')))
    #    elif mode in ['dnase_count','dnase_ave']:
    #        anno_list.append(SamHandler(pysam.Samfile(args.anno[nn])))
    #    elif mode in ['gc', 'gc_count']:
    #        anno_list.append(pysam.Fastafile(args.anno[nn]))
    #    elif mode in ['size']:
    #        anno_list.append(args.anno[nn])
    #    else:
    #        print("Unknown mode: %s" % (mode), file = sys.stderr)
    #        exit(1)
    # create center region
    for nn in tqdm(range(len(region_list))):
        region = region_list[nn]
        region.region_init(exclude_list, genome_size)
    # annotate 
    print("Annotation Region", file = sys.stderr)
    N = len(region_list)
    with Pool(args.thread) as p:
        region_list = list(tqdm(p.imap(single_region_run, region_list, 100), total = N))
    print("Annotation bed done", file = sys.stderr)
    # write to output file
    if args.output == args.bed: # re-annotate existing file
        fout = open(args.output+'.tmp', 'w')
    else:
        fout = open(args.output, 'w')
    print(region_list[0].header(), file = fout)
    if args.resort:
        region_list = sorted(region_list, key = lambda bed: (bed.chrom, bed.start))
    for region in region_list:
        try:
            print(region.write(), file = fout)
        except:
            print("[report] skip %s:%d-%d" % (region.chrom, region.start, region.end), file = sys.stderr)
            continue
    fout.close()
    # program finished 
    if args.output == args.bed:
        os.system("mv %s %s" % (args.output + '.tmp', args.output)) 

if __name__=="__main__":
    Main()
