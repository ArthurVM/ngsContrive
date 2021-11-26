import sys
import os
import subprocess
import pandas as pd

def readCSV(csv):
    df = pd.read_csv(csv, header=0)

    fout = open(f"./{os.path.splitext(os.path.basename(csv))[0]}.var", "w")
    for index, row in df.iterrows():
        varline = f"SP {row['final_annotation.Reference']} {row['final_annotation.Position']} {row['final_annotation.AlternativeNucleotide'].upper()}"
        print(varline, file=fout)

        runline = f"python3 ../ngsContrive.py -v \"{varline}\" -p {row['variant']}_ ../WD/ASM19595v2.fna /home/amorris/BioInf/art_bin_MountRainier"
        print(runline)

        subprocess.run(runline, shell=True)

    fout.close()

def main(csv):
    readCSV(csv)

if __name__=="__main__":
    if len(sys.argv) != 2:
        print("USAGE: parseResCatalogue.py <csv>")
    main(sys.argv[1])
