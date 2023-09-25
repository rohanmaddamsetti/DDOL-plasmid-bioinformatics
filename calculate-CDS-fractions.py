#!/usr/bin/env python

"""
 calculate-CDS-fractions.py by Rohan Maddamsetti.

This script answers a basic question that Lingchong asked:
for a typical plasmid or bacterial chromosome, what percentage is genuinely encoding proteins?

Usage: python calculate-CDS-fractions.py > ../results/CDS-fractions.csv
"""

import os
from Bio import SeqIO


gbk_annotation_path = "../results/gbk-annotation/"
gbk_files = [x for x in os.listdir(gbk_annotation_path) if x.endswith(".gbff")]

## print the header of the file
print("Annotation_Accession,SeqID,SeqLength,CDSLength,CDSFraction")
for gbk in gbk_files:
    AnnotationAccession = gbk.split("_genomic.gbff")[0]
    gbk_path = os.path.join(gbk_annotation_path, gbk)
    records = SeqIO.parse(gbk_path, 'genbank')
    
    ## Iterate through each SeqRecord in the GenBank file
    for record in records:
        total_sequence_length = len(record.seq)
        cds_length = 0

        ## Iterate through features and calculate CDS length
        for feature in record.features:
            if feature.type == 'CDS':
                cds_length += len(feature.location)
                
        ## Calculate the percentage of sequence covered by CDS
        CDS_fraction = cds_length / total_sequence_length
        print(",".join([AnnotationAccession, record.id, str(total_sequence_length), str(cds_length), str(CDS_fraction)]))
        
