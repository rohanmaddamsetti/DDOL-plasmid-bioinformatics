#!/usr/bin/env python

""" get_N2_symbiosis_plasmid_metadata.py by Rohan Maddamsetti. """

import csv

# Read the list of plasmids that contain Nif/Fix, Nod, or  Nol/Nop proteins.
with open('../results/plasmids_with_N2_symbiosis_genes.txt', 'r') as id_file:
    id_list = set(line.strip().replace("NC_", "").replace("NZ_", "") for line in id_file)

# Create an empty list to store the filtered rows
filtered_rows = []

# Read the chromosome and plasmid metadata file and filter rows based on the plasmid IDs
with open('../results/chromosome-plasmid-table.csv', 'r') as csv_file:
    csv_reader = csv.DictReader(csv_file)
    for row in csv_reader:
        if row['NCBI_Nucleotide_Accession'] in id_list:
            filtered_rows.append(row)

# Write the filtered rows to a new CSV file
with open('../results/plasmids_with_N2_symbiosis_metadata.csv', 'w', newline='') as output_file:
    fieldnames = ['Organism', 'Strain', 'NCBI_Nucleotide_Accession', 'SequenceType', 'Annotation_Accession']
    csv_writer = csv.DictWriter(output_file, fieldnames=fieldnames)
    csv_writer.writeheader()
    csv_writer.writerows(filtered_rows)
