## make-DDOL-figures.R by Rohan Maddamsetti.
library(tidyverse)
library(cowplot)

## ABSOLUTELY CRITICAL TODO!!! transform y-axes in plots as needed, and be mindful of ggplot pitfalls!
## ABSOLUTELY CRITICAL TODO!!!!! HANDLE ZEROS APPROPRIATELY!!!!

## Figure 6ABCD: nitrogen-fixing bacteria are enriched with metabolic genes on plasmids.

## annotate the genomes and plasmids.
gbk.annotation <- read.csv("../results/computationally-annotated-gbk-annotation-table.csv") %>%
    mutate(Annotation = replace_na(Annotation,"Unannotated")) %>%
    ## collapse Annotations into a smaller number of categories as follows:
    ## Marine, Freshwater --> Water
    ## Sediment, Soil, Terrestrial --> Earth
    ## Plants, Agriculture --> Plants
    ## Animals --> Animals
    mutate(Annotation = replace(Annotation, Annotation == "Marine", "Water")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Freshwater", "Water")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Sediment", "Earth")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Soil", "Earth")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Terrestrial", "Earth")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Plants", "Plants")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Agriculture", "Plants")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "Animals", "Animals")) %>%
    mutate(Annotation = replace(Annotation, Annotation == "blank", "Unannotated"))

## get the metadata for chromosomes and plasmids.
replicon.annotation.data <- read.csv("../results/chromosome-plasmid-table.csv") %>%
    rename(SeqType = SequenceType) %>%
    left_join(gbk.annotation) %>%
    ## Annotate the genera.
    mutate(Genus = stringr::word(Organism, 1))


plasmid.length.data <- read.csv("../results/replicon-lengths-and-protein-counts.csv") %>%
    filter(SeqType == "plasmid")

## IMPORTANT: This is a tsv file, because "," is a meaningful character in chemical names!   
plasmid.proteins.in.KEGG.metabolism <- read.table("../results/plasmid-proteins-in-KEGG-metabolism.tsv", header = TRUE)

metabolic.genes.per.plasmid <- plasmid.proteins.in.KEGG.metabolism %>%
    group_by(SeqID, SeqType) %>%
    summarize(metabolic_protein_count = n()) %>%
    mutate(metabolic_protein_count = replace_na(metabolic_protein_count, 0))


plasmid.annotation.data <- replicon.annotation.data %>%
    filter(SeqType == "plasmid")

## make the dataframe for the plot of metabolic genes on plasmids.
metabolic.gene.plot.data <- plasmid.length.data %>%
    left_join(metabolic.genes.per.plasmid) %>%
    ## make the dataframe compatible with plasmid.annotation.data
    mutate(NCBI_Nucleotide_Accession = str_remove(SeqID, "N(C|Z)_")) %>%
    ## and join.
    left_join(plasmid.annotation.data) %>%
    ## remove unannotated bacteria.
    filter(Annotation != "NA") %>%
    filter(Annotation != "Unannotated") %>%
    mutate(Annotation = factor(Annotation)) %>%
    ## Annotate Nitrogen-fixing bacteria.
    mutate(NitrogenFixer = ifelse(str_detect(Genus, regex("Rhizo|Azo", ignore_case=TRUE)), 'Nitrogen-fixer', "Non-nitrogen-fixer"))


## Figure 6ABC: plasmids with lots of metabolic proteins come from Plants and Earth.
Fig6A <- metabolic.gene.plot.data %>%
    ggplot(aes(x = metabolic_protein_count, fill = Annotation)) +
    geom_histogram(bins=100) +
    geom_vline(xintercept = 100, color="gray", linetype="dotted") +
    theme_classic() +
    scale_fill_discrete(
        drop=FALSE,
        limits = levels(metabolic.gene.plot.data$Annotation)) +
    xlim(0,500) +
    xlab("") +
    ggtitle("All plasmids with metabolic genes") +
    theme(legend.position="top")

## get the legend and remove from Fig6A.
Fig6ABC_legend <- get_legend(Fig6A)
Fig6A <- Fig6A + guides(fill = "none")

Fig6B <- metabolic.gene.plot.data %>%
    filter(metabolic_protein_count > 100) %>%
    ggplot(aes(x = metabolic_protein_count, fill = Annotation)) +
    geom_histogram(bins=100) +
    geom_vline(xintercept = 100, color="gray", linetype="dotted") +
    theme_classic() +
    scale_fill_discrete(
        drop=FALSE,
        limits = levels(metabolic.gene.plot.data$Annotation)) +
    xlim(0,500) +
    ylim(0,50) +
    xlab("") +
    ggtitle("Plasmids with > 100 metabolic genes") +
    guides(fill = "none")


Fig6C <- metabolic.gene.plot.data %>%
    filter(NitrogenFixer == "Nitrogen-fixer") %>%
    ggplot(aes(x = metabolic_protein_count, fill = Annotation)) +
    geom_histogram(bins=100) +
    geom_vline(xintercept = 100, color="gray", linetype="dotted") +
    theme_classic() +
    scale_fill_discrete(
        drop=FALSE,
        limits = levels(metabolic.gene.plot.data$Annotation)) +
    xlim(0,500) +
    ylim(0,50) +
    xlab("Number of metabolic proteins on the plasmid") +
    ggtitle("Plasmids with metabolic genes in nitrogen-fixing bacteria") +
    guides(fill = "none")
##    ggtitle("Plasmids with > 100 metabolic proteins in nitrogen-fixing bacteria") +


## metabolic genes on plasmids are enriched in nitrogen-fixing bacteria:
## graphical and statistical analysis.

N2.fixer.metabolic.gene.comparison.plot.data <- metabolic.gene.plot.data %>%
    group_by(Annotation_Accession, NitrogenFixer) %>%
    summarize(total_plasmid_metabolic_proteins = sum(metabolic_protein_count))

Fig6D <- N2.fixer.metabolic.gene.comparison.plot.data %>%
    ggplot(aes(x = total_plasmid_metabolic_proteins, y = NitrogenFixer, alpha=0.5)) +
    geom_jitter(width=0,alpha=0.2, size=0.1) +
    geom_boxplot(size=0.1, outlier.shape = NA) +
    theme_classic() +
    guides(alpha = "none") +
    xlab("Numbers of plasmid-borne metabolic genes per genome") +
    ylab("") +
    ggtitle("Plasmids in nitrogen-fixing bacteria are enriched\nwith metabolic genes")


## Save the plot
Fig6ABCD <- plot_grid(Fig6A, Fig6B, Fig6C, Fig6ABC_legend, Fig6D,
                      labels=c('A', 'B', 'C', '', 'D'), ncol=1, rel_heights=c(1,1,1,0.5,1.25))
ggsave("../results/Fig6ABCD.pdf", Fig6ABCD, width=6, height=8)


 N2.fixer.comp.data <- N2.fixer.metabolic.gene.comparison.plot.data %>%
    filter(NitrogenFixer=="Nitrogen-fixer")
non.N2.fixer.comp.data <- N2.fixer.metabolic.gene.comparison.plot.data %>%
    filter(NitrogenFixer != "Nitrogen-fixer")

## compare the means of nitrogen-fixing and non-nitrogen-fixing bacteria.
mean.N2.fixer.plasmid.metabolic.genes <- mean(N2.fixer.comp.data$total_plasmid_metabolic_proteins, na.rm=TRUE)
print(mean.N2.fixer.plasmid.metabolic.genes) ## 253.33 plasmid metabolic genes in nitrogen-fixers

mean.non.N2.fixer.plasmid.metabolic.genes <- mean(non.N2.fixer.comp.data$total_plasmid_metabolic_proteins,na.rm=TRUE)
print(mean.non.N2.fixer.plasmid.metabolic.genes) ## 22.97 plasmid metabolic genes in non-nitrogen fixers

## statistics are super significant. 
wilcox.test(
    x=N2.fixer.comp.data$total_plasmid_metabolic_proteins,
    y=non.N2.fixer.comp.data$total_plasmid_metabolic_proteins,
    alternative="greater")
            
  
