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

        if len(vl_split) != 4:
            print("\nWARNING: Weird variant line. Variant lines must be formatted as:\n\tINDEL/SP chrID 1base_index seq\n")
            exit()

        self.type = vl_split[0]
        self.chr = vl_split[1]
        self.index = int(vl_split[2])
        self.seq = vl_split[3]

def readVars(vars):
    """ read in the variant file
    """
    varbox = []

    if os.path.isfile(vars):
        with open(vars, "r") as fin:
            for line in fin.readlines():
                varbox.append(Var(line))
    else:
        varbox.append(Var(vars))

    return sorted(varbox, key=lambda x: x.index, reverse=True)

def modifySeq(fasta, vars, prefix):
    """ introduces variants defined in a variant file into sequences found in the fasta file.
    """
    seqdict = {rec.id : rec.seq for rec in SeqIO.parse(fasta, "fasta")}

    variants=0
    for var in readVars(vars):

        print(var.type, var.chr, var.index, var.seq)

        ## check index is in range
        if var.index >= len(seqdict[var.chr]):
            print(f"Index {var.index} out of range in {var.chr}.")
            continue

        newseq = seqdict[var.chr]
        real_index = var.index-1

        if var.type == "INDEL":
            variants+=1
            if var.seq.startswith("-"):
                ## deal with deletions
                newseq = newseq[:real_index] + newseq[real_index+len(var.seq):]
            else:
                ## deal with insertions
                newseq = newseq[:real_index] + var.seq + newseq[real_index:]

        elif var.type == "SP":
            ## deal with sequence polymorphisms
            variants+=1
            newseq = newseq[:real_index] + var.seq + newseq[real_index+len(var.seq):]

        seqdict[var.chr] = newseq

    fvar_path = f"./{os.path.splitext(os.path.basename(fasta))[0]}.{prefix}.fasta"
    fout = open(fvar_path, "w")
    for id, seq in seqdict.items():
        fout.write(f">{id}\n{seq}\n")

    print(f"FASTA IN : {fasta}\nVARIANTS INTRODUCED : {variants}\nFASTA OUT : {fvar_path}\n")

    # for id, seq in seqdict.items():
    #     print(id, "\t", seq)

    return fvar_path

def runART(fvar_path, artDir, prefix, doc):
    """ Run ART on the simulated variant file
    """

    artExec = os.path.join(artDir, "art_illumina")
    runline = f"{artExec} -ss HS25 -i {fvar_path} -p -l 150 -f {doc} -m 200 -s 10 -o {prefix}"

    subprocess.run(runline, shell=True, check=True, capture_output=True)

def parseArgs(argv):
    """ simple argument parser
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('fasta', action='store', help='A sequence file in FASTA format to introduce variants into and simulate reads from.')
    parser.add_argument('artDir', action='store', help='Directory containing ART executables.')

    parser.add_argument('-v', '--vars', action='store', default="", help='Either a file containing variants, or a single variant line, to introduce into the given sequence.')
    parser.add_argument('-p', '--prefix', action='store', default="varSim", help='Prefix for output files. Default=varSim.')
    parser.add_argument('-d', '--doc', action='store', default=40, help='Mean depth of coverage for simualted reads.')

    args = parser.parse_args(argv)
    return args

def main(argv):
    """ Main function
    """
    args = parseArgs(argv)

    if args.vars != "":
        fvar_path = modifySeq(args.fasta, args.vars, args.prefix)
    else:
        fvar_path = args.fasta

    runART(fvar_path, args.artDir, args.prefix, args.doc)

if __name__=="__main__":
    main(sys.argv)
