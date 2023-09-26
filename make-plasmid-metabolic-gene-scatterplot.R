## make-plasmid-metabolic-gene-scatterplot.R by Rohan Maddamsetti.
library(tidyverse)


ggplotRegression <- function(dat, xvar, yvar){
    ## code from:
    ## https://community.rstudio.com/t/annotate-ggplot2-with-regression-equation-and-r-squared/6112/7  
  fml <- paste(yvar, "~", xvar)

  fit <- lm(fml, dat)
  
  ggplot(fit$model, aes(x = .data[[xvar]], y = .data[[yvar]])) + 
      geom_point() +
      theme_classic() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}


##########################################
## data analysis starts here.

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
    left_join(gbk.annotation)

##  get output from calculate-CDS-fractions.py.
## We want to answer a basic question that Lingchong asked:
## for a typical plasmid or bacterial chromosome, what percentage is genuinely encoding proteins?
CDS.fraction.data <- read.csv("../results/CDS-fractions.csv") %>%
    ## make the dataframe compatible with replicon.annotation.data,
    mutate(NCBI_Nucleotide_Accession = str_remove(SeqID, "N(C|Z)_")) %>%
    ## and join.
    left_join(replicon.annotation.data)

####################################
## analysis of CDS fraction.

## TODO: make a binned average, so that each point is one bin of replicon size.
## Main figure, all the points together.
## supplementary figure: same figure, separated by Annotation category.

## remove NA values. haven't checked causes for this yet.
clean.CDS.fraction.data <- CDS.fraction.data %>%
    filter(!is.na(Annotation)) %>%
    filter(!is.na(SeqType))

CDS.fraction.plot1 <- ggplot(
    data = clean.CDS.fraction.data,
    aes(x = log10(SeqLength), y = log10(CDSLength), color = Annotation)) +
    geom_point(size=0.05,alpha=0.5) +
    facet_wrap(.~SeqType) +
    theme_classic()
## save the plot.
ggsave("../results/CDS-fraction-plot1.pdf", CDS.fraction.plot1, width=9)

CDS.fraction.plot2 <- ggplot(
    data = clean.CDS.fraction.data,
    aes(x = CDSFraction, fill = Annotation)) +
    geom_histogram(bins=100) +
    facet_grid(Annotation~SeqType) +
    theme_classic()
## save the plot.
ggsave("../results/CDS-fraction-plot2.pdf", CDS.fraction.plot2, width=9)


CDS.fraction.plot3 <- ggplot(
    data = clean.CDS.fraction.data,
    aes(x = log10(SeqLength), y = CDSFraction, shape = SeqType, color = SeqType)) +
    geom_point(size=0.01,alpha=0.5) +
    #facet_grid(Annotation~SeqType) +
    theme_classic()
## save the plot.
ggsave("../results/CDS-fraction-plot3.pdf", CDS.fraction.plot3, width=9)

CDS.fraction.plot4 <- ggplot(
    data = clean.CDS.fraction.data,
    aes(x = log10(SeqLength), y = CDSFraction, shape = SeqType, color = SeqType)) +
    geom_point(size=0.01,alpha=0.5) +
    facet_wrap(.~Annotation) +
    theme_classic()
## save the plot.
ggsave("../results/CDS-fraction-plot4.pdf", CDS.fraction.plot4, width=9)


####################################
## data structures for analysis of metabolic genes on plasmids.

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

## make the dataframe for the exploratory data analysis.
metabolic.gene.scatterplot.data <- plasmid.length.data %>%
    left_join(metabolic.genes.per.plasmid) %>%
    ## make the dataframe compatible with plasmid.annotation.data
    mutate(NCBI_Nucleotide_Accession = str_remove(SeqID, "N(C|Z)_")) %>%
    ## and join.
    left_join(plasmid.annotation.data)

## annotate big.plasmids.
big.plasmid.protein.threshold <- 750
big.plasmid.data <- metabolic.gene.scatterplot.data %>%
    filter(protein_count > big.plasmid.protein.threshold)

metabolic.gene.scatterplot.data <- metabolic.gene.scatterplot.data %>%
    mutate(big_plasmids = ifelse(SeqID %in% big.plasmid.data$SeqID, TRUE, FALSE))

metabolic.gene.log.scatterplot <- ggplot(
    data = metabolic.gene.scatterplot.data,
    aes(x = log2(protein_count), y = log2(metabolic_protein_count), color = big_plasmids)) +
    geom_point(size=0.2, alpha=0.5) +
    theme_classic()
## save the plot.
ggsave("../results/plasmid-metabolic-gene-log-scatterplot.pdf", metabolic.gene.log.scatterplot)

metabolic.gene.scatterplot <- ggplot(
    data = metabolic.gene.scatterplot.data,
    aes(x = protein_count, y = metabolic_protein_count, color = big_plasmids)) +
    geom_point(size=0.2, alpha=0.5) +
    theme_classic()
## save the plot.
ggsave("../results/plasmid-metabolic-gene-scatterplot.pdf", metabolic.gene.scatterplot)


metabolic.gene.log.scatterplot2 <- ggplotRegression(
    metabolic.gene.scatterplot.data,
    "log2(protein_count)", "log2(metabolic_protein_count)")
## save the plot.
ggsave("../results/plasmid-metabolic-gene-log-scatterplot2.pdf", metabolic.gene.log.scatterplot2)


metabolic.gene.log.scatterplot3 <- metabolic.gene.scatterplot.data %>%
    filter(Annotation != "blank") %>%
    filter(Annotation != "NA") %>%
    filter(Annotation != "Unannotated") %>%
    ggplot(
    aes(x = log2(protein_count), y = log2(metabolic_protein_count), color = Annotation)) +
    geom_point(size=0.2, alpha=0.5) +
    theme_classic()
## save the plot.
ggsave("../results/plasmid-metabolic-gene-log-scatterplot3.pdf", metabolic.gene.log.scatterplot3)

metabolic.gene.scatterplot3 <- metabolic.gene.scatterplot.data %>%
    filter(Annotation != "blank") %>%
    filter(Annotation != "NA") %>%
    filter(Annotation != "Unannotated") %>%
    ggplot(
    aes(x = protein_count, y = metabolic_protein_count, color = Annotation)) +
    geom_point(size=0.2, alpha=0.5) +
    theme_classic()
## save the plot.
ggsave("../results/plasmid-metabolic-gene-scatterplot3.pdf", metabolic.gene.scatterplot3)


metabolic.gene.log.scatterplot4 <- metabolic.gene.scatterplot.data %>%
    filter(Annotation != "blank") %>%
    filter(Annotation != "NA") %>%
    filter(Annotation != "Unannotated") %>%
    ggplot(
    aes(x = log2(protein_count), y = log2(metabolic_protein_count), color = Annotation)) +
    geom_point(size=0.2, alpha=0.5) +
    theme_classic() +
    facet_wrap(.~Annotation)
## save the plot.
ggsave("../results/plasmid-metabolic-gene-log-scatterplot4.pdf", metabolic.gene.log.scatterplot4)

metabolic.gene.scatterplot4 <- metabolic.gene.scatterplot.data %>%
    filter(Annotation != "blank") %>%
    filter(Annotation != "NA") %>%
    filter(Annotation != "Unannotated") %>%
    ggplot(
    aes(x = protein_count, y = metabolic_protein_count, color = Annotation)) +
    geom_point(size=0.2, alpha=0.5) +
    theme_classic() +
    facet_wrap(.~Annotation)
## save the plot.
ggsave("../results/plasmid-metabolic-gene-scatterplot4.pdf", metabolic.gene.scatterplot4)


## Super interesting. the big plasmids basically all come from nitrogen-fixing bacteria and plant pathogens!
big.plasmid.data
write.csv(x=big.plasmid.data, file="../results/big-plasmids-threshold750proteins.csv")


## calculate slope for the correlation.
metabolic.gene.log.scatterplot3 <- ggplotRegression(
    metabolic.gene.scatterplot.data,
    "protein_count", "metabolic_protein_count")
## save the plot.
ggsave("../results/plasmid-metabolic-gene-log-scatterplot3.pdf", metabolic.gene.log.scatterplot3)

## TODO:
## calculate the slope of the regression across the different ecological categories.
## look at this distribution of slope parameters, and ask whether the slope
## use a binned average, so that each range of the data gives equal contribution,
## so that this calculation does not overly favor small plasmids with few metabolic genes.

## ALSO: repeat this analysis with MGE-associated genes. examine the proportion of MGE-associated genes
## found on these plasmids across different ecological categories.
    


