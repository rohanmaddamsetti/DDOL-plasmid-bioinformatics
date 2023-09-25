#!/usr/bin/env python

"""
make-NN-plasmid-data-for-LC.py by Rohan Maddamsetti.

Usage: python make-NN-plasmid-data-for-LC.py > ../results/DDOL_plasmid_tokens.txt

Take each plasmid, and concatenate the sequence to itself 5 times.
Then, chop the 5x concatenated sequence into windows of length 512 x 3 = 1536bp.
Then, translate each window into amino acids, to generate a token that is 512aa long.

We concatenate each sequence 5 times, so that every protein on the plasmid is
contained in at least one token.

Our goal is to make a latent space that represents all known plasmids.

Another idea from LC: do the same thing for chromosomes as well?
Can we use a common embedding for plasmids and chromosomes?
If we condition on a given sequence from the user, can we predict the genomic context,
or do other interesting tasks with a downstream deep neural network?

Output specification from LC:

A *.csv file with the following columns:
NCBIPlasmidID,QuotedFASTAHeader,StartCoordinate,AAString.

The last column is a 512aa long 'peptide' that is a token generated from the plasmid.

"""

import os
from Bio.Seq import Seq


DNACHUNKLEN = 512 * 3 ## 1536 == 512 * 3.
TIMESTOROLLAROUND = 5

def get_NCBIPlasmidID(fastafilename):
    return os.path.basename(fastafilename).split(".fasta")[0]


def concatenate_the_plasmid(plasmidseq, n=TIMESTOROLLAROUND):
    return n*plasmidseq


def split_string_into_indexed_chunks(input_string, chunk_size=DNACHUNKLEN):
    original_string_length = int(len(input_string) / TIMESTOROLLAROUND)
    return [(i % original_string_length, input_string[i:i+chunk_size]) for i in range(0, len(input_string), chunk_size)]


def translate_DNA_token(DNA_token):
    return str(Seq(DNA_token).translate())

def generate_NN_plasmid_data(plasmids_fasta_dir):
    ## print the header.
    print("NCBIPlasmidID,QuotedFASTAHeader,StartCoordinate,AAString")

    ## now go through all plasmids in complete genomes as of August 2023.
    fasta_filelist = [x for x in os.listdir(plasmids_fasta_dir) if x.endswith(".fasta")]
    for fasta_file in fasta_filelist:
        NCBIPlasmidID = get_NCBIPlasmidID(fasta_file)
        fasta_path = os.path.join(plasmids_fasta_dir, fasta_file)
        with open(fasta_path) as plasmid_fasta_fh:
            header = plasmid_fasta_fh.readline().strip()
            plasmidseq = plasmid_fasta_fh.readline().strip()
            quoted_header = "\"" + header + "\""
            concatenated_plasmidseq = concatenate_the_plasmid(plasmidseq, n=5)
            ## remove any dangling tokens with length != 512.
            untranslated_token_tuples = [(i,tok) for i, tok in split_string_into_indexed_chunks(concatenated_plasmidseq) if len(tok) == DNACHUNKLEN]
            token_tuples = [(token_index, translate_DNA_token(DNA_token)) for token_index, DNA_token in untranslated_token_tuples]
            for token_index, token in token_tuples:
                row = ','.join([NCBIPlasmidID, quoted_header, str(token_index), token])
                print(row)
    return


def main():
    plasmids_fasta_dir = "../results/plasmid-FASTA-references/"
    generate_NN_plasmid_data(plasmids_fasta_dir)
    return


if __name__ == "__main__":
    main()
