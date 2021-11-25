import sys
import os
import argparse
from Bio import SeqIO

class Var(object):
    """docstring for Var."""

    def __init__(self, varline):
        super(Var, self).__init__()
        vl_split = varline.strip("\n").split(" ")
        self.type = vl_split[0]
        self.chr = vl_split[1]
        self.index = int(vl_split[2])
        self.seq = vl_split[3]

def readVarFile(variantFile):
    """ read in the variant file
    """
    varbox = []

    with open(variantFile, "r") as fin:
        for line in fin.readlines():
            varbox.append(Var(line))

    return sorted(varbox, key=lambda x: x.index)

def modifySeq(fasta, variantFile):
    """ introduces variants defined in a variant file into sequences found in the fasta file.
    """
    seqdict = {rec.id : rec.seq for rec in SeqIO.parse(fasta, "fasta")}

    offset=0
    for var in readVarFile(variantFile):

        print(var.type, var.chr, var.index, var.seq)

        ## check index is in range
        if var.index >= len(seqdict[var.chr]):
            print(f"Index {var.index} out of range in {var.chr}.")
            continue

        newseq = str(seqdict[var.chr])

        if var.type == "INDEL":
            if var.seq.startswith("-"):
                ## deal with deletions
                offset-=len(var.seq)
                newseq = newseq[:var.index-1+offset] + newseq[:var.index-1+offset+len(var.seq)]
            else:
                ## deal with insertions
                offset+=len(var.seq)
                newseq = newseq[:var.index-1+offset] + var.seq + newseq[:var.index-1+offset]

        elif var.type == "SNP":
            newseq[var.index+offset-1] == var.seq

        seqdict[var.chr] = newseq

    for id, seq in seqdict.items():
        print(id, seq)

def parseArgs(argv):
    """ simple argument parser
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('fasta', action='store', help='A sequence file in FASTA format to introduce variants into and simulate reads from.')
    parser.add_argument('variantFile', action='store', help='A file containing variants to introduce into the given sequence.')
    parser.add_argument('artDir', action='store', help='Directory containing ART executables.')

    parser.add_argument('-p', '--prefix', action='store', default="varSim", help='Prefix for output files. Default=varSim.')

    args = parser.parse_args(argv)
    return args

def main(argv):
    """ Main function
    """
    args = parseArgs(argv)
    modifySeq(args.fasta, args.variantFile)

if __name__=="__main__":
    main(sys.argv)
