#!/usr/bin/env python

"""
DDOL.py by Maggie Wilson and Rohan Maddamsetti.

This script calculates the entropy of the distribution of genes
in the Nif/Fix Nitrogen fixation pathway, the Nod nodulation gene pathway,
and the Nol/Nop effector protein pathway,
defined by the genes in Figure 3 of Weisberg et al. (2022) in mBio:
'Pangenome Evolution Reconciles Robustness and Instability of Rhizobial Symbiosis'
and the nolXWBTUV cultivar-specificity locus based on Deakin and Broughton (2009)
"Symbiotic use of pathogenic strategies: rhizobial protein secretion systems" in Nat. Rev. Micro.

"""

import csv
import re
import math

###########
## Functions

def count_plasmids_containing_pathway_genes(pathway_plasmid_proteins_file):
    seqid_list = []
    with open(pathway_plasmid_proteins_file, "r") as pathway_plasmids_fh:
        for line in pathway_plasmids_fh:
            seqid_pattern = r'>SeqID=([^\|]+)'
            seqid_match = re.search(seqid_pattern, line)
            if seqid_match:
                ## get the plasmid SeqID, and remove NZ_ or NC_ prefixes.
                seqid = seqid_match.group(1).replace("NZ_","").replace("NC_","")
                ## we are implicitly counting the number of genes in the given pathway
                ## found on a particular plasmid-- this list needs to have duplicate
                ## entries for each time a gene in the pathway is found on the plasmid!
                seqid_list.append(seqid)
    return (seqid_list)


def get_genome_accessions(infile="../results/genomes_with_N2_symbiosis_plasmid_genes.csv"):
    annotation_acc_list = []
    with open(infile, "r") as infile_fh:
        csv_candidate_reader = csv.DictReader(infile_fh, delimiter=',')
        for row in csv_candidate_reader:
            annotation_acc_list.append(row['Annotation_Accession'])
    ## return a sorted list of unique genome IDs.
    unique_annotation_list = sorted(list(set(annotation_acc_list)))
    print(f"analyzing {len(unique_annotation_list)} candidate DDOL genomes")
    return (unique_annotation_list)


def make_plasmid_to_genome_dict(candidate_genome_list,
                                chromosome_plasmid_table="../results/chromosome-plasmid-table.csv"):
    plasmid_to_genome_dict = {}
    with open(chromosome_plasmid_table, "r") as chr_plasmid_table_file:
        csv_plasmid_reader = csv.DictReader(chr_plasmid_table_file, delimiter=',')
        for row in csv_plasmid_reader:
            genome_ann_acc = row.get('Annotation_Accession','')
            sequence_type = row.get('SequenceType','')
            plasmid_id = row.get('NCBI_Nucleotide_Accession')
            ## only consider plasmids in the candidate DDOL genomes of interest.
            if genome_ann_acc in candidate_genome_list and sequence_type == 'plasmid':
                plasmid_to_genome_dict[plasmid_id] = genome_ann_acc
    return (plasmid_to_genome_dict)

    

def initialize_genome_to_plasmid_Pathway_dict(candidate_genome_list,
                                              chromosome_plasmid_table="../results/chromosome-plasmid-table.csv"):
    genome_to_plasmid_Pathway_dict = {}
    with open(chromosome_plasmid_table, "r") as chr_plasmid_table_file:
        csv_plasmid_reader = csv.DictReader(chr_plasmid_table_file, delimiter=',')
        for row in csv_plasmid_reader:
            genome_ann_acc = row.get('Annotation_Accession','')
            sequence_type = row.get('SequenceType','')
            plasmid_id = row.get('NCBI_Nucleotide_Accession')
            ## only consider genomes in the candidate DDOL genomes of interest.
            if genome_ann_acc in candidate_genome_list and sequence_type == 'plasmid':      
                genome_to_plasmid_Pathway_dict.setdefault(genome_ann_acc, {})[plasmid_id] = 0
    return (genome_to_plasmid_Pathway_dict)


def make_genome_to_plasmid_Pathway_dict(unique_annotation_list,
                                        plasmid_to_genome_dict,
                                        plasmid_pathway_count_list):
    genome_to_plasmid_pathway_dict = initialize_genome_to_plasmid_Pathway_dict(unique_annotation_list)

    ## hack to match plasmid IDs when the prefixes don't match.
    potential_plasmids_in_plasmid_to_genome_dict = plasmid_to_genome_dict.keys()
    
    for plasmid in plasmid_pathway_count_list:
        ## match potential plasmids in plasmid_to_genome_dict.
        plasmid_matches1 = [x for x in potential_plasmids_in_plasmid_to_genome_dict if re.search(plasmid, x)]
        if not len(plasmid_matches1): ## skip if the plasmid is not found
            continue
        matching_plasmid1 = plasmid_matches1[0]
        my_genome = plasmid_to_genome_dict[matching_plasmid1]
        ## only consider genomes in the set of candidate DDOL genomes of interest.
        if my_genome in genome_to_plasmid_pathway_dict:
            ## match potential plasmids to get the ID in the dictionary.
            potential_plasmids_in_cur_genome = genome_to_plasmid_pathway_dict[my_genome].keys()
            plasmid_matches2 = [x for x in potential_plasmids_in_cur_genome if re.search(plasmid, x)]
            matching_plasmid2 = plasmid_matches2[0] ## assume there's only one match.
            ## now update the closest matching plasmid in the given genome in the data structure.
            genome_to_plasmid_pathway_dict[my_genome][matching_plasmid2] += 1
    return (genome_to_plasmid_pathway_dict)
    

def calculate_plasmid_probabilities(genome_to_plasmid_dict):
    prob_dict = {}
    for genome_ann_acc, plasmid_count_dict in genome_to_plasmid_dict.items():
        total_count = sum(plasmid_count_dict.values())
        for plasmid_id, count in plasmid_count_dict.items():
            if total_count > 0:
                probability = float(count)/float(total_count)
            else:
                probability = 0
            prob_dict.setdefault(genome_ann_acc, {})[plasmid_id] = probability
    return prob_dict


def calculate_entropy(prob_dict):
    entropy_dict = {}
    for genome_ann_acc, probabilities in prob_dict.items():
        entropy = 0
        for plasmid_id, probability in probabilities.items():
            if probability > 0:
                entropy += (-probability) * math.log2(probability)
            entropy_dict[genome_ann_acc] = entropy
    return entropy_dict


def write_entropy_dict(combined_entropy_dict, outfile="../results/entropy_data.csv"):
    with open(outfile, 'w', newline='') as entropy_file:
        fieldnames = ['Annotation_Accession','Entropy','PathwayType']
        writer = csv.DictWriter(entropy_file, fieldnames=fieldnames)
        writer.writeheader()
        for source, entropy_dict in combined_entropy_dict.items():
            for genome_ann_acc, entropy in entropy_dict.items():
                ##if isinstance(entropy, int) and entropy == 0 and not isinstance(entropy, float):
                writer.writerow({
                    'Annotation_Accession': genome_ann_acc,
                    'Entropy': entropy,
                    'PathwayType': source
                })


def print_pathway_gene_count_dict(genome_to_plasmid_pathway_dict, pathwaytype):
    assert pathwaytype in ["Nif/Fix", "Nod", "All Symbiosis Pathways"]
    for genome, plasmid_dict in genome_to_plasmid_pathway_dict.items():
        total_pathway_gene_count = sum(plasmid_dict.values())
        ## print out interesting genomes
        if total_pathway_gene_count > 0:
            print(f"Genome: {genome}: Plasmid distribution of {pathwaytype}: {plasmid_dict}")

                
def write_pathway_gene_count_dict_to_csv(genome_to_plasmid_pathway_dict, outfile, pathwaytype):
    assert pathwaytype in ["Nif/Fix", "Nod", "All Symbiosis Pathways"]
    with open(outfile, "w") as outfh:
        print("Annotation_Accession,Plasmid,PathwayType,PathwayGeneCount", file=outfh)
        for genome, plasmid_dict in genome_to_plasmid_pathway_dict.items():
            for plasmid, pathway_gene_count in plasmid_dict.items():
                print(f"{genome},{plasmid},{pathwaytype},{pathway_gene_count}", file=outfh)
    return


def main():

    ## make lists of the plasmids containing genes in each of the three pathways of interest.
    NifFix_plasmid_file = "../results/nif_fix_plasmid_proteins.faa"
    plasmid_nif_count_list = count_plasmids_containing_pathway_genes(NifFix_plasmid_file)

    Nod_plasmid_file = "../results/nod_plasmid_proteins.faa"
    plasmid_nod_count_list = count_plasmids_containing_pathway_genes(Nod_plasmid_file)

    combined_pathway_count_list = plasmid_nif_count_list + plasmid_nod_count_list

    ## get a sorted list of the genomes (Annotation Accessions).
    candidate_genome_list = get_genome_accessions()

    ## make a dictionary of plasmid ID to genome ID (Annotation Accession)
    plasmid_to_genome_dict = make_plasmid_to_genome_dict(candidate_genome_list)

    ## For each pathway in [Nod, Nif/Fix]:
    #	Make a dict of dicts with the following structure:
    #	genome_to_plasmid_NodPathway_dict = { genomeID: {plasmidID1: 0, plasmidID2 : 0   … , plasmidID7 : 0} }

    ## make the big data structures.
    genome_to_plasmid_NifFix_Pathway_dict = make_genome_to_plasmid_Pathway_dict(candidate_genome_list,
                                                                                plasmid_to_genome_dict, plasmid_nif_count_list)

    genome_to_plasmid_Nod_Pathway_dict = make_genome_to_plasmid_Pathway_dict(candidate_genome_list,
                                                                             plasmid_to_genome_dict, plasmid_nod_count_list)

    genome_to_plasmid_all_Pathways_dict = make_genome_to_plasmid_Pathway_dict(candidate_genome_list,
                                                                              plasmid_to_genome_dict, combined_pathway_count_list)

    ## write these dictionaries into csv files.
    write_pathway_gene_count_dict_to_csv(genome_to_plasmid_NifFix_Pathway_dict, "../results/NifFixPlasmidCount.csv", "Nif/Fix")
    write_pathway_gene_count_dict_to_csv(genome_to_plasmid_Nod_Pathway_dict, "../results/NodPlasmidCount.csv", "Nod")
    write_pathway_gene_count_dict_to_csv(genome_to_plasmid_all_Pathways_dict, "../results/AllPathwaysPlasmidCount.csv", "All Symbiosis Pathways")
    
    ## print out interesting genomes
    print_pathway_gene_count_dict(genome_to_plasmid_NifFix_Pathway_dict, "Nif/Fix")
    print()
    print_pathway_gene_count_dict(genome_to_plasmid_Nod_Pathway_dict, "Nod")
    print()
    print_pathway_gene_count_dict(genome_to_plasmid_all_Pathways_dict, "All Symbiosis Pathways")
            
    ## For each of the pathways, iterate through the 
    ## 	genome_to_plasmid_Pathway_dict data structure, and calculate 
    ## 	the entropy of the distribution of genes across its plasmids. 
    ## 	Save the result in a dictionary,  genome_to_pathway_entropy_dict
    NifFix_prob = calculate_plasmid_probabilities(genome_to_plasmid_NifFix_Pathway_dict)
    Nod_prob = calculate_plasmid_probabilities(genome_to_plasmid_Nod_Pathway_dict)
    combined_pathways_prob = calculate_plasmid_probabilities(genome_to_plasmid_all_Pathways_dict)

    NifFix_entropy = calculate_entropy(NifFix_prob)
    Nod_entropy = calculate_entropy(Nod_prob)
    combined_pathways_entropy = calculate_entropy(combined_pathways_prob)

    combined_entropy_dict = {'Nif/Fix':NifFix_entropy, 'Nod':Nod_entropy, "All Symbiosis Pathways": combined_pathways_entropy}

    write_entropy_dict(combined_entropy_dict)
    return


## run the main function.
if __name__ == "__main__":
    main()
    
