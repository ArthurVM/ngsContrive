import os, sys
import math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def subsample(fasta, perc):
    """ takes a single sequence fasta and extracts the middle portion of it totalling the given percentage of the whole genome
    """
    outfile = f"{os.path.basename(os.path.splitext(fasta)[0])}.{perc}SS.fa"
    with open(outfile, "w") as fout:
        for rec in SeqIO.parse(fasta, "fasta"):
            chunk = (perc/100)*len(rec.seq)
            l = math.floor(chunk/2)
            r = math.ceil(chunk/2)
            m = math.floor(len(rec.seq)/2)
            subsample = rec.seq[m-l:m+r]
            record = SeqRecord(
                Seq(subsample),
                id=f"{rec.id}_{m-l}:{m+r}",
                name=rec.id,
                description=f"| {m-l}:{m+r}")

            SeqIO.write(record, fout, "fasta")

def main(argv):
    fasta = argv[1]
    perc = argv[2]
    subsample(fasta, int(perc))

if __name__=="__main__":
    if len(sys.argv) != 3:
        print("USAGE : subsampleFASTA.py <fasta> <sample_percentage>")
        sys.exit()
    main(sys.argv)
