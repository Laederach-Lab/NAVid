import argparse, os, subprocess, pickle
from pathlib import Path
from Bio import SeqIO

link = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz"

# arguments
parser = argparse.ArgumentParser(prog="parse-hg38.py", description=f"Parse reference genome sequence data. If FASTA file is not provided, it will be downloaded from {link}")
parser.add_argument("-f", "--fasta", metavar="FASTA", type=Path, required=False,
                    help=f"path to GRCh38/hg38 sequence file; can be zipped")
args = parser.parse_args()

ROOT = Path(__file__).parent
OUT_DIR = f"{ROOT}/nagnag/data"
os.makedirs(f"{OUT_DIR}/seq", exist_ok=True)

# Download
if args.fasta == None:
    print(f"Downloading {link}")
    hg38_fa_gz = Path(f"{ROOT}/hg38.fa.gz")
    subprocess.run(["wget", "-O", hg38_fa_gz, link], check=True)
else:
    hg38_fa_gz = Path(args.fasta)
    if not hg38_fa_gz.is_file():
        parser.error("Invalid reference genome FASTA!")

# Decompress
if hg38_fa_gz.suffix == ".gz":
    hg38_fa = Path(f"{ROOT}/hg38.fa")
    print(f"Unzipping {hg38_fa_gz} to {hg38_fa}")
    subprocess.run(["gunzip", hg38_fa_gz, "-c"], stdout=open(hg38_fa, "w"), check=True)
else:
    hg38_fa = hg38_fa_gz

# Parse
print(f"Reading {hg38_fa}")
hg38_dict = {record.id : record.seq for record in SeqIO.parse(hg38_fa, "fasta")}

canon_fmt = [f"chr{c}" for c in list(map(str, range(1, 23))) + ["X", "Y", "M"]]

# Canonical Chromosomes (dictionary)
hg38 = hg38_dict.copy()
keys = list(hg38)
[hg38.pop(c) for c in keys if c not in canon_fmt]
pickle.dump(hg38, open(f"{OUT_DIR}/hg38", "wb"))

# By Chromosome
hg38 = hg38_dict.copy()

# canonical chromosomes
for chrom in canon_fmt:
    pickle.dump(hg38.pop(chrom), open(f"{OUT_DIR}/seq/{chrom}", "wb"))

# non-canonical chromosomes
pickle.dump(hg38, open(f"{OUT_DIR}/seq/non-canonical", "wb"))
