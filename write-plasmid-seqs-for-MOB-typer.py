#!/usr/bin/env python

"""
write-plasmid-seqs-for-MOB-typer.py by Rohan Maddamsetti.

This script writes out plasmids as FASTA sequences into ../results/plasmid-FASTA-seqs
to make inputs for MOB-typer (a program in MOB-SUITE).
"""

import os
from tqdm import tqdm
from Bio import SeqIO

def writeout_all_plasmids_as_fasta(refgenomes_dir, FASTA_plasmid_dir):
    gbk_filelist = [x for x in os.listdir(refgenomes_dir) if x.endswith("gbff")]
    for gbk_file in tqdm(gbk_filelist):
        gbk_path = os.path.join(refgenomes_dir, gbk_file)
        genome_id = gbk_file.split(".gbff")[0]
        with open(gbk_path,'r') as gbk_fh:
            SeqID = None
            SeqType = None
            for record in SeqIO.parse(gbk_fh, "genbank"):
                SeqID = record.id
                if "plasmid" in record.description:
                    SeqType = "plasmid"
                else:
                    continue
                header = ">" + "|".join(["SeqID="+SeqID,"SeqType="+SeqType,"description="+record.description])
                ## open the file for writing.
                outfilename = genome_id + "_" + SeqID + ".fasta"
                outfile = os.path.join(FASTA_plasmid_dir, outfilename)
                with open(outfile, "w") as outfh:
                    outfh.write(header + "\n")
                    outfh.write(str(record.seq) + "\n")
    return
    

def main():
    refgenomes_dir = "../results/gbk-annotation/"
    FASTA_plasmid_dir = "../results/plasmid-FASTA-references"
    """ make plasmids as input to MOB-typer."""
    writeout_all_plasmids_as_fasta(refgenomes_dir, FASTA_plasmid_dir)


## run the script.
main()
