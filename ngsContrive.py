import sys
import os
import argparse
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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

    variants=0
    offset=0
    for var in readVarFile(variantFile):

        # print(var.type, var.chr, var.index, var.seq)

        ## check index is in range
        if var.index >= len(seqdict[var.chr]):
            print(f"Index {var.index} out of range in {var.chr}.")
            continue

        newseq = seqdict[var.chr]
        real_index = var.index+offset

        if var.type == "INDEL":
            variants+=1
            if var.seq.startswith("-"):
                ## deal with deletions
                offset-=len(var.seq)
                newseq = newseq[:real_index] + newseq[real_index+len(var.seq):]
            else:
                ## deal with insertions
                offset+=len(var.seq)
                newseq = newseq[:real_index] + var.seq + newseq[real_index:]

        elif var.type == "SNP":
            variants+=1
            newseq = newseq[:real_index-1] + var.seq + newseq[real_index:]

        seqdict[var.chr] = newseq

    fvar_path = f"./{os.path.splitext(os.path.basename(fasta))[0]}.var.fasta"
    fout = open(fvar_path, "w")
    for id, seq in seqdict.items():
        fout.write(f">{id}\n{seq}\n")

    print(f"FASTA IN : {fasta}\nVARIANTS INTRODUCED : {variants}\nFASTA OUT : {fvar_path}\n")

    return fvar_path

def runART(fvar_path, artDir, prefix):
    """ Run ART on the simulated variant file
    """

    artExec = os.path.join(artDir, "art_illumina")
    runline = f"{artExec} -ss HS25 -i {fvar_path} -p -l 150 -f 40 -m 200 -s 10 -o {prefix}"

    subprocess.run(runline, shell=True, check=True, capture_output=True)

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
    fvar_path = modifySeq(args.fasta, args.variantFile)
    runART(fvar_path, args.artDir, args.prefix)

if __name__=="__main__":
    main(sys.argv)
