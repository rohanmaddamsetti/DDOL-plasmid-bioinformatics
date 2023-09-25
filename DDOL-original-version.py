#!/usr/bin/env python

"""
DDOL.py by Maggie Wilson and Rohan Maddamsetti.

This script calculates the entropy of the distribution of genes
in the Nif/Fix Nitrogen fixation pathway, the Nod nodulation gene pathway,
and the T3SS/other nodulation effector protein pathway,
defined by the genes in Figure 3 of Weisberg et al. (2022) in mBio:
'Pangenome Evolution Reconciles Robustness and Instability of Rhizobial Symbiosis'

"""

import csv
import re
import math


## 1) For each pathway in [Nod, Nif/Fix, T3SS/symbiosis effectors]:
#	Initialize a dict of dicts with the following structure:
#	genome_to_plasmid_NodPathway_dict = { genomeID: {plasmidID1: 0, plasmidID2 : 0   … , plasmidID7 : 0} }
annotation_acc_list = []
plasmid_id = []

genome_to_plasmid_Pathway_dict = {}

## Make sure these path names are correct when running on your local machine.
candidate_csv_file = "../results/candidate_DDOL_genomes_with_N2_symbiosis_or_T3SS_genes.csv"
chromosome_plasmid_table = "../results/chromosome-plasmid-table.csv"
NifFix_plasmid_file = "../results/nif_fix_plasmid_proteins.faa"
Nod_plasmid_file = "../results/nod_plasmid_proteins.faa"
T3SS_plasmid_file = "../results/nodulation_and_T3SS_effector_plasmid_proteins.faa"

with open(candidate_csv_file, "r") as candidate_file:
    csv_candidate_reader = csv.DictReader(candidate_file, delimiter=',')
    for row in csv_candidate_reader:
        annotation_acc_list.append(row['Annotation_Accession'])
        annotation_set = set(annotation_acc_list)
        unique_annotation_list = list(annotation_set)
        
with open(chromosome_plasmid_table, "r") as chr_plasmid_table_file:
    csv_plasmid_reader = csv.DictReader(chr_plasmid_table_file, delimiter=',')
    for row in csv_plasmid_reader:
        plasmid_ann_acc = row.get('Annotation_Accession','')
        sequence_type = row.get('SequenceType','')
        plasmid_id = row.get('NCBI_Nucleotide_Accession')
        
        if plasmid_ann_acc in unique_annotation_list and sequence_type == 'plasmid':
            plasmid_counts = 0
            
            genome_to_plasmid_Pathway_dict.setdefault(plasmid_ann_acc, {})[plasmid_id] = plasmid_counts
print(genome_to_plasmid_Pathway_dict)

seqid_nif_list =[]
with open(NifFix_plasmid_file, "r") as NifFix_file:
    for line in NifFix_file:
        seqid_pattern = r'>SeqID=([^\|]+)'
        seqid_nif_match = re.search(seqid_pattern, line)
        if seqid_nif_match:
            seqid_nif = seqid_nif_match.group(1)
            seqid_nif = seqid_nif[3:]
            seqid_nif_list.append(seqid_nif)
print(seqid_nif_list)

seqid_nod_list =[]
with open(Nod_plasmid_file, "r") as Nod_file:
    for line in Nod_file:
        seqid_pattern = r'>SeqID=([^\|]+)'
        seqid_nod_match = re.search(seqid_pattern, line)
        if seqid_nod_match:
            seqid_nod = seqid_nod_match.group(1)
            seqid_nod = seqid_nod[3:]
            seqid_nod_list.append(seqid_nod)
print(seqid_nod_list)

seqid_T3SS_list =[]
with open(T3SS_plasmid_file, "r") as T3SS_file:
    for line in T3SS_file:
        seqid_pattern = r'>SeqID=([^\|]+)'
        seqid_T3SS_match = re.search(seqid_pattern, line)
        if seqid_T3SS_match:
            seqid_T3SS = seqid_T3SS_match.group(1)
            seqid_T3SS = seqid_T3SS[3:]
            seqid_T3SS_list.append(seqid_T3SS)
print(seqid_T3SS_list)

genome_to_plasmid_NifFix_Pathway_dict = {}
genome_to_plasmid_Nod_Pathway_dict = {}
genome_to_plasmid_T3SS_Pathway_dict = {}

for plasmid_ann_acc, plasmid_dict in genome_to_plasmid_Pathway_dict.items():
    for plasmid_id, plasmid_counts in plasmid_dict.items():
        count_NifFix = seqid_nif_list.count(plasmid_id)
        count_Nod = seqid_nod_list.count(plasmid_id)
        count_T3SS = seqid_T3SS_list.count(plasmid_id)
        
        genome_to_plasmid_Pathway_dict[plasmid_ann_acc][plasmid_id]= count_NifFix + count_Nod + count_T3SS
        if plasmid_id in seqid_nif_list:
            genome_to_plasmid_NifFix_Pathway_dict.setdefault(plasmid_ann_acc, {})[plasmid_id]=count_NifFix
        elif plasmid_id not in seqid_nif_list:
            genome_to_plasmid_NifFix_Pathway_dict.setdefault(plasmid_ann_acc, {})[plasmid_id]=0
            
        if plasmid_id in seqid_nod_list:
            genome_to_plasmid_Nod_Pathway_dict.setdefault(plasmid_ann_acc, {})[plasmid_id]=count_Nod
        elif plasmid_id not in seqid_nod_list:
            genome_to_plasmid_Nod_Pathway_dict.setdefault(plasmid_ann_acc, {})[plasmid_id]=0
            
        if plasmid_id in seqid_T3SS_list:
            genome_to_plasmid_T3SS_Pathway_dict.setdefault(plasmid_ann_acc, {})[plasmid_id]=count_T3SS
        elif plasmid_id not in seqid_T3SS_list:
            genome_to_plasmid_T3SS_Pathway_dict.setdefault(plasmid_ann_acc, {})[plasmid_id]=0


print(genome_to_plasmid_Pathway_dict)
print(genome_to_plasmid_NifFix_Pathway_dict)
print(genome_to_plasmid_Nod_Pathway_dict)
print(genome_to_plasmid_T3SS_Pathway_dict)


## 4) For each of the pathways, iterate through the 
# 	genome_to_plasmid_Pathway_dict data structure, and calculate 
# 	the entropy of the distribution of genes across its plasmids. 
# 	Save the result in a dictionary,  genome_to_pathway_entropy_dict

NifFix_prob = {}
Nod_prob = {}
T3SS_prob = {}

def calculate_plasmid_probabilities(genome_to_plasmid_dict):
    prob_dict = {}
    for plasmid_ann_acc, plasmid_counts in genome_to_plasmid_dict.items():
        total_count = sum(plasmid_counts.values())
        for plasmid_id, count in plasmid_counts.items():
            if total_count > 0:
                probability = float(count)/float(total_count)
            else:
                probability = 0
            prob_dict.setdefault(plasmid_ann_acc, {})[plasmid_id] = probability
    return prob_dict

NifFix_prob = calculate_plasmid_probabilities(genome_to_plasmid_NifFix_Pathway_dict)
Nod_prob = calculate_plasmid_probabilities(genome_to_plasmid_Nod_Pathway_dict)
T3SS_prob = calculate_plasmid_probabilities(genome_to_plasmid_T3SS_Pathway_dict) 
        
        
print(NifFix_prob)
print(Nod_prob)
print(T3SS_prob)



NifFix_entropy = {}
Nod_entropy = {}
T3SS_entropy = {}

def calculate_entropy(prob_dict):
    entropy_dict = {}
    for plasmid_ann_acc, probabilities in prob_dict.items():
        entropy = 0
        for plasmid_id, probability in probabilities.items():
            if probability > 0:
                entropy += (-probability) * math.log2(probability)
            entropy_dict[plasmid_ann_acc] = entropy
    return entropy_dict

NifFix_entropy = calculate_entropy(NifFix_prob)
Nod_entropy = calculate_entropy(Nod_prob)
T3SS_entropy = calculate_entropy(T3SS_prob)

print(NifFix_entropy)
print(Nod_entropy)
print(T3SS_entropy)

combined_entropy_dict = {'NifFix':NifFix_entropy, 'Nod':Nod_entropy, 'Nodulation and T3SS Effectors':T3SS_entropy}
combined_prob_dict = {'NifFix':NifFix_prob, 'Nod':Nod_prob, 'Nodulation and T3SS Effectors':T3SS_prob}

with open('../results/entropy_data.csv', 'w', newline='') as entropy_file:
    fieldnames = ['Plasmid Annotation Accession','Entropy','Pathway Type']
    writer = csv.DictWriter(entropy_file, fieldnames=fieldnames)
    writer.writeheader()
    for source, entropy_dict in combined_entropy_dict.items():
        for plasmid_ann_acc, entropy in entropy_dict.items():
            if isinstance(entropy, int) and entropy == 0 and not isinstance(entropy, float):
                writer.writerow({
                'Plasmid Annotation Accession': plasmid_ann_acc,
                'Entropy': entropy,
                'Pathway Type': source
                })


    
#####################################################################################
#Extra un-needed code
#####################################################################################


################ This section not completely needed
NifFix_sum = {}
Nod_sum = {}
T3SS_sum = {}

def sum_plasmid_counts_per_genome(genome_to_plasmid_dict):
    genome_counts = {}
    for plasmid_ann_acc, plasmid_counts in genome_to_plasmid_dict.items():
        total_count = sum(plasmid_counts.values())
        #print(total_count)
        genome_counts[plasmid_ann_acc] = total_count
        print(genome_counts)
    return genome_counts

NifFix_sum = sum_plasmid_counts_per_genome(genome_to_plasmid_NifFix_Pathway_dict)
Nod_sum = sum_plasmid_counts_per_genome(genome_to_plasmid_Nod_Pathway_dict)
T3SS_sum = sum_plasmid_counts_per_genome(genome_to_plasmid_T3SS_Pathway_dict)

print(NifFix_sum)
print(Nod_sum)
print(T3SS_sum)

#######################


## 2) Create a dictionary mapping Plasmid IDs 
# 	to the GenomeID of the genome containing that plasmid:
#	plasmid_to_genome_dict = {plasmidID1 : genomeID1, plasmidID2 : genomeID1, plasmidID3 : genomeID2, … }
## 3) For each of the pathways, iterate through the corresponding file in 
# 	[nodulation_and_T3SS_effector_plasmid_proteins.faa, nod_plasmid_proteins.faa, 
# 	nif_fix_plasmid_proteins.faa], and update the corresponding 
# 	genome_to_plasmid_Pathway_dict data structure.

plasmid_to_genome_NifFix_Pathway_dict = {}
plasmid_to_genome_Nod_Pathway_dict = {}
plasmid_to_genome_T3SS_Pathway_dict = {}

for plasmid_ann_acc, plasmid_dict in genome_to_plasmid_NifFix_Pathway_dict.items():
    for plasmid_id, value in plasmid_dict.items():
        plasmid_to_genome_NifFix_Pathway_dict[value] = plasmid_ann_acc
        
    
for plasmid_ann_acc, plasmid_dict in genome_to_plasmid_Nod_Pathway_dict.items():
    for plasmid_id, value in plasmid_dict.items():
        plasmid_to_genome_Nod_Pathway_dict[value] = plasmid_ann_acc
        
        
for plasmid_ann_acc, plasmid_dict in genome_to_plasmid_T3SS_Pathway_dict.items():
    for plasmid_id, value in plasmid_dict.items():  
        plasmid_to_genome_T3SS_Pathway_dict[value] = plasmid_ann_acc 
        
print(plasmid_to_genome_NifFix_Pathway_dict)
print(plasmid_to_genome_Nod_Pathway_dict)
print(plasmid_to_genome_T3SS_Pathway_dict)

