## find_N2_symbiosis_plasmid_genomes.R by Rohan Maddamsetti
library(tidyverse)

DDOL.genome.initial.candidates <- read.csv("../results/plasmids_with_N2_symbiosis_metadata.csv")

DDOL.candidates <- DDOL.genome.initial.candidates %>%
    filter(SequenceType == "plasmid") %>%
    group_by(Organism, Strain, Annotation_Accession) %>%
    summarize(NumberOfPlasmidsWith_N2_symbiosis_genes = n()) %>%
    filter(NumberOfPlasmidsWith_N2_symbiosis_genes >= 1) %>%
    arrange(desc(NumberOfPlasmidsWith_N2_symbiosis_genes))

DDOL.genome.final.candidates <- DDOL.genome.initial.candidates %>%
    filter(Annotation_Accession %in% DDOL.candidates$Annotation_Accession)

write.csv(DDOL.genome.final.candidates, "../results/genomes_with_N2_symbiosis_plasmid_genes.csv",
          quote=F, row.names=FALSE)
