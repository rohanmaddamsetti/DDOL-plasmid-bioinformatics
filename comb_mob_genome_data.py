#!/usr/bin/env python

""" comb_mob_genome_data.py by Hye-in Son.
combines the mob typer results into a single .csv file
"""

import glob

# Create a new '.txt' file for writing
with open('../results/combined_mob_typer_results.txt','w',newline='') as txt_file):
	first_file = True # Flag to track the first file because I want the header to show only once

	# Iterate through each '.txt' file in the data directory (plasmid-MOB-typer-results-rerun/)
	for path in glob.glob("../plasmid-MOB-typer-results-rerun/*results.txt"):
		with open("../plasmid-MOB-typer-results-rerun/"+path, 'r') as txt_input_file:
			# Read the content of the input text file
			txt_content = txt_input_file.read()

			# Write the heading only in the first file
			if first_file:
				txt_file.write(txt_content)
				first_file = False
			else:
				# Skip the heading in subsequent files
				lines = txt_content.split('\n')[1:]
				txt_file.write('\n'.join(lines))
