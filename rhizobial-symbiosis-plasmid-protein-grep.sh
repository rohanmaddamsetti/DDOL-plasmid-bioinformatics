#!/usr/bin/env bash

## rhizobial-symbiosis-plasmid-protein-grep.sh by Rohan Maddamsetti

## This script uses the gene/protein names of symbiosis genes
## involved in nitrogen fixation in Rhizobia bacteria that
## are found in SynICE MGEs, reported in Weisberg et al. (2022)
## in MBio: "Pangenome Evolution Reconciles Robustness and Instability of Rhizobial Symbiosis"

## The protein names used in these regexes are found in the columns of the matrix
## in Figure 3.

cat ../results/plasmid-protein-db.faa | grep -A 1 -E "Fix[ABCRIXUS]|Nif[ABDEHKNQSWXZ]" | grep -v '^--$' > ../results/nif_fix_plasmid_proteins.faa
cat ../results/plasmid-protein-db.faa | grep -A 1 -E "Nod[ABCSUIJDFMNYZVW]|Noe[DIL]|Nol[AKMNOYZ]" | grep -v '^--$' > ../results/nod_plasmid_proteins.faa
cat ../results/plasmid-protein-db.faa | grep -A 1 -E "Nol[BUV]|T3SS|type_3_secretion_system|type_III_secretion_system|Nop[ABDELMPTXYZ]" | grep -v '^--$' > ../results/nodulation_and_T3SS_effector_plasmid_proteins.faa

## Now, get all unique plasmid SeqIDs matched in these three files.
cat ../results/nif_fix_plasmid_proteins.faa ../results/nod_plasmid_proteins.faa ../results/nodulation_and_T3SS_effector_plasmid_proteins.faa | sed -nr 's/.*SeqID=([^|]*)\|.*/\1/p' | sort | uniq > ../results/plasmids_with_N2_symbiosis_or_T3SS_genes.txt

## Now get the genomic metadata for these plasmids.
python get_N2_symbiosis_and_T3SS_plasmid_metadata.py

## Now, get all genomes that have multiple plasmids that matched.
Rscript find_N2_symbiosis_and_T3SS_DDOL_candidates.R

## output is saved in: candidate_DDOL_genomes_with_N2_symbiosis_or_T3SS_genes.csv