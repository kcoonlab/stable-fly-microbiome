set.seed(123456)

source("./Helper_Functions.R")

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

variable1 = as.character(get_variable(ps, "sample.type"))
variable2 = as.character(get_variable(ps, "location"))
variable3 = as.character(get_variable(ps, "sampling.date"))
sample_data(ps)$NewPastedVar <- mapply(paste0, variable1, variable2, variable3, collapse = "_")
ps.Merged <- merge_samples(ps, "NewPastedVar")
metadata <- data.frame(sample_data(ps.Merged))
orderAbdDF <- MakeAbundanceDF(physeq = ps.Merged,
                               taxRank = "Order",
                               abundanceFilter = 0.02)
#write.csv(orderAbdDF, "orderAbdDF.csv", row.names=FALSE)
        
orderAbdDF_formatted <- read.csv('orderAbdDF_formatted.csv',header=TRUE)
metadata <- metadata %>%
  rownames_to_column(var = "Sample")
sampleLevels <- as.character(metadata$Sample)
sampleLabels <- metadata %>%
  dplyr::select(Sample)
sampleLabelsVec <- sampleLabels$Sample
names(sampleLabelsVec) <- sampleLabels$Sample
orderAbdDF_formatted$Sample <- factor(orderAbdDF_formatted$Sample, levels = sampleLevels)
orderAbdDF_formatted$Order <- factor(orderAbdDF_formatted$Order, levels=c(" Acetobacterales"," Caulobacterales"," Rhizobiales"," Rhodobacterales"," Sphingomonadales",
									" Aeromonadales"," Burkholderiales"," Cardiobacteriales"," Chromatiales"," Enterobacterales",
									" Pseudomonadales"," Xanthomonadales"," Bacillales"," Exiguobacterales"," Lactobacillales",
									" Staphylococcales"," Christensenellales"," Clostridiales"," Clostridia"," Eubacteriales",
									" Lachnospirales"," Monoglobales"," Oscillospirales"," Peptostreptococcales-Tissierellales"," Acholeplasmatales",
									" Erysipelotrichales"," Flavobacteriales"," Bacteroidales"," Bifidobacteriales"," Corynebacteriales",
									" Micrococcales"," Propionibacteriales"," Coriobacteriales"," Campylobacterales"," Spirochaetales"))
cbPalette <- unique(orderAbdDF_formatted$Color)
orderCompPlot <- ggplot(orderAbdDF_formatted,
                         aes_string(x = "Sample", y = "Abundance",
                                    fill = "Order")) +
  geom_bar(stat = "identity", width = 1) +
  theme_bw() +
  scale_x_discrete(labels = sampleLabelsVec) +
  scale_fill_manual(values=cbPalette)
orderCompPlot








