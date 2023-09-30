#!/usr/bin/env python

""" join_candidate_DDOL_mob.py by Hye-in Son.

Rohan's plasmids_with_N2_symbiosis_metadata.csv has 5 columns:
	- Organism
	- Strain
	- NCBI_Nucleotide_Accession
	- SequenceType
	- Annotation_Accession

My mob_typer results has a column (sample_id) that is basically a combination of Rohan's "Annotation_Accession" + "_genomic_NZ_" + "NCBI_Nucleotide_Accession".

This code is modifying Rohan's so that the new output file has a new column that matches the format of my results' "sample_id" column.
Eventually, I will join the modified file of Rohan's with mine.
"""

import csv

# Define the input and output file names
input_file = '../results/plasmids_with_N2_symbiosis_metadata.csv'
output_file = '../results/DDOL_symbiosis_plasmids.csv'

# Open the original CSV file for reading and the modified file for writing
with open(input_file, 'r', newline='') as infile, open(output_file, 'w', newline='') as outfile:
    # Create CSV reader and writer objects
    csv_reader = csv.reader(infile)
    csv_writer = csv.writer(outfile)

    # Read the header and add the new column header
    header = next(csv_reader)
    header.append('sample_id')  # Add the new column header
    csv_writer.writerow(header)  # Write the modified header to the output file

    # Process and write the remaining rows
    for row in csv_reader:
        # Extract the values from the columns you want to merge
        col1_value, col2_value = row[4], row[2]

        # Merge the values from col1 and col2 into a single string
        merged_value = col1_value  + "_genomic_NZ_" + col2_value

        # Append the merged value as a new column to the row
        row.append(merged_value)

        # Write the modified row to the output file
        csv_writer.writerow(row)

print(f'Modified file "{output_file}" has been created.')
