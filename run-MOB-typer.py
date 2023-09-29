#!/usr/bin/env python

"""
run-MOB-typer.py by Rohan Maddamsetti and Hye-in Son.

This script writes out plasmids as FASTA sequences into ../results/plasmid-FASTA-seqs
to make inputs for MOB-typer (a program in MOB-SUITE).

This script is run on DCC inside of a wrapper shell script written by Hye-in as follows:
sbatch run-MOB-typer.sh
"""

import os
from tqdm import tqdm
import subprocess


def run_MOB_typer(FASTA_plasmid_dir, MOB_typer_results_dir):
    plasmid_files = [x for x in os.listdir(FASTA_plasmid_dir) if x.endswith(".fasta")]
    for plasmidfile in tqdm(plasmid_files):
        plasmid_path = os.path.join(FASTA_plasmid_dir, plasmidfile)
        MOB_typer_outfile = plasmidfile.split(".fasta")[0] + "_mobtyper_results.txt"
        MOB_typer_outpath = os.path.join(MOB_typer_results_dir, MOB_typer_outfile)
        MOB_typer_args = ["mob_typer", "-i", plasmid_path, "-o", MOB_typer_outpath]
        print (" ".join(MOB_typer_args))
        subprocess.run(MOB_typer_args)
    return


def main():
    FASTA_plasmid_dir = "../results/plasmid-FASTA-references/"
    MOB_typer_results_dir = "../results/plasmid-MOB-typer-results/"
    run_MOB_typer(FASTA_plasmid_dir, MOB_typer_results_dir)
    return


## run the script.
main()
