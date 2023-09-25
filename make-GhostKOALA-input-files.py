#!/usr/bin/env python

"""
make-GhostKOALA-input-files.py by Rohan Maddamsetti.

This script splits ../results/plasmid-protein-db.faa into pieces that are <300MB with <500,000 sequences per file,
so that each piece can be submitted to the GhostKOALA webserver at: https://www.kegg.jp/ghostkoala/

Since I have 4,157,061 plasmid proteins in total, split into pieces with 462,000 sequences each, so that we have 9 files in total.
NOTE: each sequence also has a associated header, so there is a factor of 2 in the script to account for this.
"""

import os

MAX_SEQUENCES = 462000
MAX_LINES = 2 * MAX_SEQUENCES


def makeGhostKOALA_input_filename(cur_num):
    GhostKOALA_inputfile_dir = "../results/GhostKOALA-input-files/"
    return os.path.join(GhostKOALA_inputfile_dir, "GhostKOALA_input_" + str(cur_num) + ".fasta")


def split_file_into_chunks_by_max_lines(input_file_path, max_lines=MAX_LINES):
    with open(input_file_path, 'r') as input_file:
        cur_num = 1
        ## fancy way of removing empty lines from the list in the comprehension.
        line_buffer = [x for x in [input_file.readline() for i in range(max_lines)] if len(x)]
        while len(line_buffer): ## while the buffer is non-empty.
            output_file_path = makeGhostKOALA_input_filename(cur_num)
            with open(output_file_path, 'w') as output_file:
                output_file.writelines(line_buffer)
            cur_num += 1
            ## fancy way of removing empty lines from the list in the comprehension.
            line_buffer = [x for x in [input_file.readline() for i in range(max_lines)] if len(x)]
    return


def main():
    big_plasmid_protein_dbfile = "../results/plasmid-protein-db.faa"
    split_file_into_chunks_by_max_lines(big_plasmid_protein_dbfile)
    return

## run the script.
if __name__ == "__main__":
    main()
