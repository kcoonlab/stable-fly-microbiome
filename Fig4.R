set.seed(123456)

# Load required packages

# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")

library(phyloseq)
library(dplyr)
library(tidyverse)

# Read in phyloseq object, ensure data format is consistent across sample metadata
ps <- readRDS("ps_FieldWork2021_AJS_Final.rds")
sample_data(ps)$SampleID <- row.names(sample_data(ps))
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="7.9.21"] <- "7.9.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="7.16.21"] <- "7.16.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="7.23.21"] <- "7.23.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="7.30.21"] <- "7.30.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="8.6.21"] <- "8.6.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="8.13.21"] <- "8.13.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="8.20.21"] <- "8.20.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="8.27.21"] <- "8.27.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="9.10.21"] <- "9.10.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="7.8.21"] <- "7.8.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="7.15.21"] <- "7.15.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="7.22.21"] <- "7.22.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="7.29.21"] <- "7.29.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="8.5.21"] <- "8.5.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="8.12.21"] <- "8.12.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="8.19.21"] <- "8.19.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="8.26.21"] <- "8.26.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="9.2.21"] <- "9.2.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="9.9.21"] <- "9.9.2021"
sample_data(ps)$sampling.date[sample_data(ps)$sampling.date=="9.15.21"] <- "9.15.2021"
ps

# Drop samples with less than 100 total reads
ps <- subset_samples(ps, sample_sums(ps)>100) 

# Drop the ASVs that have zeros across all samples
ps <- subset_taxa(ps,taxa_sums(ps)>0)

# Drop ASVs assigned to "Archaea", "Mitochondria" or "Chloroplast"
ps <- ps %>% subset_taxa( Domain!= "Archaea" & Family!= " Mitochondria" | is.na(Family) )
ps <- ps %>% subset_taxa( tax_table(ps)[,"Order"]!=" Chloroplast" | is.na(tax_table(ps)[,"Order"]) )

# Collect sample data
SamDat = data.frame(sample_data(ps))
SamDat$SampleID = row.names(SamDat)

ps.norm = transform_sample_counts(ps, function(x) x / sum(x) )
ps = ps.norm

#### PCoA to compare communities in different sample types from each facility ####
## Arlington samples first

# Create PCoA ordination
ps.arlington <- subset_samples(ps, location=="Arlington")
ps.ordination.PCoA = ordinate(ps.arlington, method="PCoA", distance="bray")

# plot ordination
p = plot_ordination(ps.arlington, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p = p + theme_bw() + geom_point(size =1.5)
p = p + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34", #Endo
                                      "#19468C", #Ecto
                                      "#B75463"),#Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p = p + theme(legend.text = element_text(size = 10))
p$layers = p$layers[-1] #to remove the larger point coded in original plot
p

## Now DCC samples

# Create PCoA ordination
ps.dcc <- subset_samples(ps, location=="DCC")
ps.ordination.PCoA = ordinate(ps.dcc, method="PCoA", distance="bray")

# plot ordination
p = plot_ordination(ps.dcc, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p = p + theme_bw() + geom_point(size =1.5)
p = p + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34", #Endo
                                      "#19468C", #Ecto
                                      "#B75463"),#Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p = p + theme(legend.text = element_text(size = 10))
p$layers = p$layers[-1] #to remove the larger point coded in original plot
p
