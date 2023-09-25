#!/usr/bin/env bash

## rhizobial-symbiosis-plasmid-protein-grep.sh by Rohan Maddamsetti

## This script uses the gene/protein names of symbiosis genes
## involved in nitrogen fixation in Rhizobia bacteria that
## are found in SynICE MGEs, reported in Weisberg et al. (2022)
## in MBio: "Pangenome Evolution Reconciles Robustness and Instability of Rhizobial Symbiosis"
## I also added the nolXWBTUV cultivar-specificity locus based on Deakin and Broughton (2009)
## "Symbiotic use of pathogenic strategies: rhizobial protein secretion systems" in Nat. Rev. Micro.

## Usage: ./rhizobial-symbiosis-plasmid-protein-grep.sh

## The protein names used in these regexes are found in the columns of the matrix
## in Figure 3. I also added nolXWBTUV from Deakin et al. (2009).

cat ../results/plasmid-protein-db.faa | grep -A 1 -E "Fix[ABCRIXUS]|Nif[ABDEHKNQSWXZ]" | grep -v '^--$' > ../results/nif_fix_plasmid_proteins.faa
cat ../results/plasmid-protein-db.faa | grep -A 1 -E "Nod[ABCSUIJDFMNYZVW]|Noe[DIL]|Nol[AKMNOYZ]" | grep -v '^--$' > ../results/nod_plasmid_proteins.faa
cat ../results/plasmid-protein-db.faa | grep -A 1 -E "Nol[BUVXWT]|Nop[ABDELMPTXYZ]" | grep -v '^--$' > ../results/nol_nop_plasmid_proteins.faa

## Now, get all unique plasmid SeqIDs matched in these three files.
cat ../results/nif_fix_plasmid_proteins.faa ../results/nod_plasmid_proteins.faa ../results/nol_nop_plasmid_proteins.faa | sed -nr 's/.*SeqID=([^|]*)\|.*/\1/p' | sort | uniq > ../results/plasmids_with_N2_symbiosis_genes.txt

## Now get the genomic metadata for these plasmids.
python get_N2_symbiosis_plasmid_metadata.py

## Now, get all genomes that have one or more plasmids that matched.
Rscript find_N2_symbiosis_plasmid_genomes.R
## output is saved in: genomes_with_N2_symbiosis_plasmid_genes.csv

## make input files for DDOL_entropy.R.
python DDOL.py

## Now I can run the analysis in DDOL_entropy.R.
