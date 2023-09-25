#!/usr/bin/env python

"""
get_unique_KEGG_IDs.py by Rohan Maddamsetti.

This script gets all unique KEGG IDs in ../results/filtered_successful_GhostKOALA_concatenated_ko_results.tsv,
and writes a file with a placeholder first column, and the unique KEGG IDs in the second column,
per the format required by KEGG Mapper Reconstruct: https://www.genome.jp/kegg/mapper/reconstruct.html

Usage: python get_unique_KEGG_IDs.py > ../results/unique_plasmid_KEGG_IDs.tsv

Then, upload unique_plasmid_KEGG_IDs.tsv to KEGG Mapper Reconstruct:
https://www.genome.jp/kegg/mapper/reconstruct.html
"""


unique_KEGG_set = set()
with open("../results/filtered_successful_GhostKOALA_concatenated_ko_results.tsv") as infh:
    for line in infh:
        line = line.strip() ## removing leading/lagging whitespace.
        header, KEGG_ID = line.split("\t")
        unique_KEGG_set.add(KEGG_ID)

unique_KEGG_list = sorted(list(unique_KEGG_set))

for i,KEGG_ID in enumerate(unique_KEGG_list):
    ## we need a placeholder in the first column to make KEGG Mapper Reconstruct happy.
    print(f"x{i}\t{KEGG_ID}")
