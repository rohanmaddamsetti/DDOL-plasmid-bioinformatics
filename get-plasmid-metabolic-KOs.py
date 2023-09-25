#!/usr/bin/env python

"""
get-plasmid-metabolic-KOs.py by Rohan Maddamsetti.

Usage: python get-plasmid-metabolic-KOs.py

I submitted ../results/unique_plasmid_KEGG_IDs.tsv to the KEGG Mapper
Reconstruct webserver, and saved a list of plasmid 2265 KO IDs that map to
KEGG Pathway 01100 "Metabolic pathways", in the file ../results/unique_plasmid_metabolic_KEGG_IDs.txt.

This script filters ../results/filtered_GhostKOALA_concatenated_ko_results.tsv
for plasmid proteins with KOs in the list in
../results/unique_plasmid_metabolic_KEGG_IDs.txt.

"""


def get_replicons_from_genomes_with_chromosomes_smaller_than_a_plasmid(replicon_length_file):
    bad_replicon_list  = []
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
                if (SeqType == "plasmid") and (replicon_length > cur_chromosome_length) and (SeqID not in bad_replicon_list):
                    bad_replicon_list.append(SeqID)
    return bad_replicon_list


def main():

    """ IMPORTANT: This is a tsv file, because "," is a meaningful character in chemical names!   """
    replicon_length_file = "../results/replicon-lengths-and-protein-counts.csv"
    """ make a list of genomes in which the chromosome is smaller than some plasmid.
    These are often misannotated genomes in which the real chromosome is switched up with a plasmid."""
    bad_replicon_list = get_replicons_from_genomes_with_chromosomes_smaller_than_a_plasmid(replicon_length_file)

    print("*******")
    print("BAD REPLICONS FROM GENOMES WITH A PLASMID APPARENTLY LARGER THAN THE CHROMOSOME:")
    for x in bad_replicon_list:
        print(x)
    print("*******")

    ## for debugging other scripts: what is the set of SeqID prefixes?
    NCBI_prefixes = set()
    
    metabolic_KEGG_IDs = set()
    with open("../results/unique_plasmid_metabolic_KEGG_IDs.txt", "r") as metabolic_KEGG_fh:
        for line in metabolic_KEGG_fh:
            metabolic_KEGG_ID = line.strip()
            metabolic_KEGG_IDs.add(metabolic_KEGG_ID)

    line_buffer = [] ## This is to save the lines to write to file.
    with open("../results/successful_GhostKOALA_concatenated_ko_results.tsv") as all_plasmid_GhostKOALA_fh:
        for line in all_plasmid_GhostKOALA_fh:
            line = line.strip()
            plasmid_gene_string, KO_ID = line.split("\t")
            SeqID, SeqType, locus_tag, product, = [ x.split("=")[-1] for x in plasmid_gene_string.split("|") ]
            NCBI_prefix = SeqID.split("_")[0]
            NCBI_prefixes.add(NCBI_prefix)

            ## skip over replicons in genomes in which there is a plasmid larger than the chromosome.
            ## In some cases this is due to misannotation of chromosomes and plasmids in RefSeq genomes.
            if SeqID in bad_replicon_list:
                continue
            if KO_ID in metabolic_KEGG_IDs:
                line_buffer.append("\t".join([SeqID, SeqType, locus_tag, product]))

    ##  write line_buffer to the output file.
    with open("../results/plasmid-proteins-in-KEGG-metabolism.tsv", "w") as outfh:
        outfh.write("SeqID\tSeqType\tlocus_tag\tproduct\n")
        for line in line_buffer:
            outfh.write(line + "\n")

    ## print the set of NCBI prefixes.
    print("THE SET OF NCBI PREFIXES:")
    print(NCBI_prefixes)
    print("*******")

    return


if __name__ == "__main__":
    main()
