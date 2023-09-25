## find_N2_symbiosis_and_T3SS_DDOL_candidates.R by Rohan Maddamsetti


library(tidyverse)

DDOL.genome.initial.candidates <- read.csv("../results/plasmids_with_N2_symbiosis_or_T3SS_genes_metadata.csv")

DDOL.candidates <- DDOL.genome.initial.candidates %>%
    filter(SequenceType == "plasmid") %>%
    group_by(Organism, Strain, Annotation_Accession) %>%
    summarize(NumberOfPlasmidsWithCandidateDDOL_N2_symbiosis_or_T3SS_genes = n()) %>%
    filter(NumberOfPlasmidsWithCandidateDDOL_N2_symbiosis_or_T3SS_genes >= 2) %>%
    arrange(desc(NumberOfPlasmidsWithCandidateDDOL_N2_symbiosis_or_T3SS_genes))

DDOL.genome.final.candidates <- DDOL.genome.initial.candidates %>%
    filter(Annotation_Accession %in% DDOL.candidates$Annotation_Accession)

write.csv(DDOL.genome.final.candidates, "../results/candidate_DDOL_genomes_with_N2_symbiosis_or_T3SS_genes.csv",
          quote=F, row.names=FALSE)
