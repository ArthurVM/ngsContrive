import sys
import os
from Bio import SeqIO

def modifySeq(fasta, variantFile):
    """ introduces variants defined in a variant file into sequences found in the fasta file.
    """


def parseArgs(argv):
    """ simple argument parser
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('fasta', action='store', help='A sequence file in FASTA format to introduce variants into and simulate reads from.')
    parser.add_argument('variantFile', action='store', help='A file containing variants to introduce into the given sequence.')
    parser.add_argument('artDir', action='store', help='Directory containing ART executables.')

    parser.add_argument('-p', 'prefix', action='store', default="varSim", help='Prefix for output files. Default=varSim.')

    args = parser.parse_args(argv)
    return args

def main(argv):
    """ Main function
    """



if __name__=="__main__":
    main(sys.argv)
