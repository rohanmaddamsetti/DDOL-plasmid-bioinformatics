#!/usr/bin/env python

"""
make-plasmid-protein-FASTA-db.py by Rohan Maddamsetti.

This script goes through all the gbk annotation files in ../results/gbk-annotation,
and makes a FASTA database of all proteins found on plasmids.

"""

import os
from tqdm import tqdm
from Bio import SeqIO


def find_genomes_with_chromosomes_smaller_than_a_plasmid(replicon_length_file):
    bad_genomes_list  = []
    cur_annotation_accession = None
    cur_chromosome_length = -1
    with open(replicon_length_file, "r") as replicon_length_fh:
        for i, line in enumerate(replicon_length_fh):
            if i == 0: continue ## skip the header
            line = line.strip() ## remove leading and lagging whitespace.
            AnnotationAccession, SeqID, SeqType, replicon_length, protein_count = line.split(",")
            replicon_length = int(replicon_length)
            protein_count = int(protein_count)
            if AnnotationAccession != cur_annotation_accession and SeqType == "chromosome":
                cur_annotation_accession = AnnotationAccession
                cur_chromosome_length = replicon_length
            else: ## AnnotationAccession == cur_annotation_accession
                if (SeqType == "plasmid") and (replicon_length > cur_chromosome_length) and (cur_annotation_accession not in bad_genomes_list):
                    bad_genomes_list.append(cur_annotation_accession)
    return bad_genomes_list


def generate_plasmid_protein_fasta_db(refgenomes_dir, bad_genomes_list):
    """
    IMPORTANT: This function only writes out proteins that are found on plasmids.
    """
    
    fasta_outfile = "../results/plasmid-protein-db.faa"
    with open(fasta_outfile, "w") as outfh:
        gbk_filelist = [x for x in os.listdir(refgenomes_dir) if x.endswith("gbff")]
        good_gbk_filelist = [x for x in gbk_filelist if x not in bad_genomes_list]
        for gbk_file in tqdm(good_gbk_filelist):
            gbk_path = os.path.join(refgenomes_dir, gbk_file)
            genome_id = gbk_file.split(".gbff")[0]

            AnnotationAccession = genome_id.split("_genomic")[0]
            """ skip genomes with a chromosome smaller than a plasmid in the genome.
            This is probably a RefSeq genome with a 'plasmid' that is actually the chromosome.
            """ 
            if AnnotationAccession in bad_genomes_list: continue
            
            with open(gbk_path,'r') as gbk_fh:
                SeqID = None
                SeqType = None
                for record in SeqIO.parse(gbk_fh, "genbank"):
                    SeqID = record.id
                    if "plasmid" in record.description:
                        SeqType = "plasmid"
                    else:
                        continue
                    for feature in record.features:
                        ## only analyze protein-coding genes.
                        if feature.type != "CDS": continue
                        ##print(feature.qualifiers)
                        locus_tag = feature.qualifiers["locus_tag"][0]
                        ## for downstream parsing, replace spaces with underscores in the product annotation field.
                        product = feature.qualifiers["product"][0].replace(" ","_")
                        ## skip over CDS that don't have an annotated translation.
                        if "translation" not in feature.qualifiers: continue
                        protein_seq = feature.qualifiers["translation"][0]
                        header = ">" + "|".join(["SeqID="+SeqID,"SeqType="+SeqType,"locus_tag="+locus_tag,"product="+product])
                        outfh.write(header + "\n")
                        outfh.write(str(protein_seq) + "\n")
    return


def main():
    refgenomes_dir = "../results/gbk-annotation/"
    replicon_length_file = "../results/replicon-lengths-and-protein-counts.csv"

    """ make a list of genomes in which the chromosome is smaller than some plasmid.
    These are often misannotated genomes in which the real chromosome is switched up with a plasmid."""
    bad_genomes_list = find_genomes_with_chromosomes_smaller_than_a_plasmid(replicon_length_file)

    print("*******")
    print("BAD GENOMES WITH A PLASMID APPARENTLY LARGER THAN THE CHROMOSOME:")
    for x in bad_genomes_list:
        print(x)
    print("*******")
    
    """ make a database of proteins found on plasmids as input to GhostKOALA to get KEGG IDs.  """
    generate_plasmid_protein_fasta_db(refgenomes_dir, bad_genomes_list)
    return


## run the script.
if __name__ == "__main__":
    main()
