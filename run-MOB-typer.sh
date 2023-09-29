#!/bin/bash
## run-MOB-typer.sh by Hye-in Son.
## Usage: sbatch run-MOB-typer.sh

#SBATCH --output=output/script-%j.out
#SBATCH --mem=32G
#SBATCH -p youlab-gpu
#SBATCH --cpus-per-task=1
#SBATCH --job-name=jobname-%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hyein.son@duke.edu

python run-MOB-typer.py
