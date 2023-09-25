#!/usr/bin/env python

""" remove-chromosomes-from-Ghost_KOALA-results.py by Rohan Maddamsetti.

This script is a temporary downstream hack to remove proteins found on the sequences in bad-replicons.txt
from successful_GhostKOALA_concatenated_ko_results.tsv.

Usage: python remove-chromosomes-from-Ghost_KOALA-results.py > ../results/filtered_successful_GhostKOALA_concatenated_ko_results.tsv

"""

bad_replicon_list = []

with open("../results/bad-replicons.txt", "r") as bad_replicon_fh:
    for line in bad_replicon_fh:
        if line.startswith("#"): continue ## skip headers
        bad_replicon = line.strip()
        bad_replicon_list.append(bad_replicon)


with open("../results/successful_GhostKOALA_concatenated_ko_results.tsv", "r") as data_fh:
    for line in data_fh:
        line = line.strip()
        SeqID = line.split('|')[0].split('=')[-1]
        this_one_bad = False
        for x in bad_replicon_list:
            if x == SeqID:
                this_one_bad = True
                break
        if not this_one_bad:
            print(line)
