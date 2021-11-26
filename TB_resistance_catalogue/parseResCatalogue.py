import sys
import os
import pandas as pd

def readCSV(csv):
    df = pd.read_csv(csv, header=0)

    fout = open(f"./{os.path.splitext(os.path.basename(csv))[0]}.var", "w")
    for index, row in df.iterrows():
        print(f"SP NC_000962.3 {row['genome_index']} {row['alt_nt'].upper()}", file=fout)
    fout.close()

def main(csv):
    readCSV(csv)

if __name__=="__main__":
    if len(sys.argv) != 2:
        print("USAGE: parseResCatalogue.py <csv>")
    main(sys.argv[1])
