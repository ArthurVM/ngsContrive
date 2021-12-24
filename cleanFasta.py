import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(fasta):
    b = ["A", "C", "T", "G", "N"]
    newrecs = []
    c= 0

    for rec in SeqIO.parse(fasta, "fasta"):
        newseq = ""
        for n in rec.seq.upper():
            if n not in b:
                c += 1
                newseq += "N"
            else:
                newseq += n
        record = SeqRecord(Seq(newseq), id=rec.id, name=rec.id, description=rec.description)
        newrecs.append(record)

    SeqIO.write(newrecs, f"{fasta}.cleaned", "fasta")
    print(f"Wrote new fasta with {c} substitutions.")

if __name__=="__main__":
    if len(sys.argv) != 2:
        print("USAGE : cleanFasta.py <fasta>")
        sys.exit()
    main(sys.argv[1])
