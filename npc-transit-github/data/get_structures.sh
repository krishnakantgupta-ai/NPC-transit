#!/bin/sh
# Downloads the PDB entries used in data/dataset_mohr2009.csv into data/pdb/ under the file names the dataset expects.
# Needs internet access (curl) and awk. Where an entry contains several copies, only chain A is kept.
set -e
cd "$(dirname "$0")"
mkdir -p pdb
get() { curl -sf -o "pdb/$2.pdb" "https://files.rcsb.org/download/$1.pdb" || { echo "download failed: $1"; exit 1; }; }
chainA() { awk '$0 !~ /^(ATOM|HETATM|ANISOU|TER)/ || substr($0,22,1)=="A"' "pdb/$1.pdb" > "pdb/$1.tmp" && mv "pdb/$1.tmp" "pdb/$1.pdb"; }
get 2JV1 INSULIN        # human insulin monomer, NMR (two chains A and B are both kept)
get 5PTI 5PTI           # aprotinin
get 1FIL PROFILIN1      # human profilin 1
get 1UBQ 1UBQ           # ubiquitin
get 2SPZ 2SPZ           # Z-domain of protein A (NMR)
get 2TRX 2TRX; chainA 2TRX               # E. coli thioredoxin (reference probe, not scored)
get 1F6S LACTALBUMIN; chainA LACTALBUMIN # bovine alpha-lactalbumin, chain A
get 1EMA 1EMA           # GFP (wild type; the assay used enhanced GFP)
get 1IXH PBP            # E. coli phosphate-binding protein with phosphate
get 1ANF 1ANF           # E. coli maltose-binding protein
echo "Structures written to data/pdb/. Check titles with: grep -m1 '^TITLE' data/pdb/*.pdb"
