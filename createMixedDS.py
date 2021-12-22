import sys
import re
import argparse
from math import ceil
from os import path
from collections import defaultdict
from Bio import SeqIO

def isFile(filename):
    """ Checks if a path is a file """

    if not path.isfile(filename):
        time = strftime("%H:%M:%S", localtime())
        print(
            f"{time} :: No file found at {filename}",
            end="\n",
            file=sys.stderr,
            flush=True,
        )
        exit(3)
    else:
        return path.abspath(path.realpath(path.expanduser(filename)))

def isDir(dirname):
    """ Checks if a path is a directory """

    if not path.isdir(dirname):
        time = strftime("%H:%M:%S", localtime())
        print(
            f"{time} :: No directory found at {dirname}",
            end="\n",
            file=sys.stderr,
            flush=True,
        )
        exit(3)
    else:
        return path.abspath(path.realpath(path.expanduser(dirname)))

def parseMixFile(mix_file):
    """ Parse the mix file
    """
    prop_dict = defaultdict(float)
    with open(mix_file, "r") as fin:
        for line in fin.readlines():
            id, prop = line.split(" ")
            prop_dict[id] = float(prop)

    sumprop = sum([prop for id, prop in prop_dict.items()])
    if round(sumprop, 4) != 1:
        print(f"The summed proportions do not equal 1 ({sumprop}). Please check the mix file.")
        sys.exit(1)

    return prop_dict

def getReads(fqpath, numreads):
    """ gets n reads from a fastq file
    """
    fqbox = []
    for i, read in enumerate(SeqIO.parse(fqpath, "fastq")):
        if i == numreads:
            break
        else:
            fqbox.append(read)

    return fqbox

def flatten(t):
    return [item for sublist in t for item in sublist]

def genDataset(prop_dict, input_dir, total_reads, out):
    forward_reads = []
    reverse_reads = []

    for id, p in prop_dict.items():
        numreads = ceil((total_reads*p)*1000)
        # print(f"Getting {numreads} reads from {id}")
        f_path = path.join(input_dir, f"{id}_1.fq")
        r_path = path.join(input_dir, f"{id}_2.fq")

        forward_reads.append(getReads(f_path, numreads))
        reverse_reads.append(getReads(r_path, numreads))

    SeqIO.write(flatten(forward_reads), f"{out}_1.fq", "fastq")
    SeqIO.write(flatten(reverse_reads), f"{out}_2.fq", "fastq")

def parseArgs(argv):
    """ simple argument parser
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('mix_file', action='store', type=isFile, help='A file containing contaminant read proportions.')
    parser.add_argument('input_dir', action='store', type=isDir, help='Input directory containing contaminant paired end read files.')
    parser.add_argument('total_reads', action='store', type=int, help='The total number of reads for dataset construction (in 10^3).')
    parser.add_argument('-o', '--out_prefix', action='store', default=sys.stdout, help='Prefix for output file.')

    args = parser.parse_args(argv)
    return args

def main(argv):
    args = parseArgs(argv)
    prop_dict = parseMixFile(args.mix_file)
    genDataset(prop_dict, args.input_dir, args.total_reads, args.out_prefix)

if __name__=="__main__":
    main(sys.argv)
