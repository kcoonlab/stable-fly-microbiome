set.seed(123456)

# Load required packages

# install.packages("devtools")
# devtools::install_github("adw96/breakaway")
library(breakaway)
# You may need to update Breakaway and then close and reopen R

# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")

# if (!requireNamespace("BiocManager", quietly = TRUE)){
#    install.packages("BiocManager")
# }
# BiocManager::install("ANCOMBC")

library(phyloseq)
library(tidyverse)
library(dplyr)
library(purrr)
library(ANCOMBC)

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

#### Generate Breakaway Richness Estimates ####

# Set OTU table
otu_data = t(otu_table(ps))
# Set sample data
meta_data = sample_data(ps)
# Flip OTU table so rownames match sample data
head(colnames(otu_data) == rownames(meta_data))

# Run Breakaway's frequency table list function
frequencytablelist = build_frequency_count_tables(otu_data)

# Run the richness estimator (breakaway) on all our samples (lists of frequency tables)
RichEsts = lapply(frequencytablelist,breakaway)

# Pull out the estimates, errors, and the model
Estimate = as.matrix(map(RichEsts, "estimate"))
Error = as.matrix(map(RichEsts, "error"))
Model = as.matrix(map(RichEsts, "model"))
df = data.frame(Estimate,Error,Model)

# Add sample ID column, estimate, and error
df$SampleID = row.names(df)
df$Estimate=as.numeric(df$Estimate)
df$Error=as.numeric(df$Error)

#Remove rows where error is equal to zero
df = df[df$Error>0.01,]

# Merge the estimates with the sample data
RichPlot3 = merge(SamDat,df,by="SampleID")
head(RichPlot3)

#### Summarize the estimates across different sample types (internal flies, external flies, manure) within each facility (Arlington, DCC) ####

# First, make a single variable that has both of those elements
RichPlot3$Comp = paste(RichPlot3$location,RichPlot3$sample.type)
RichPlot3$Comp = as.factor(RichPlot3$Comp)
head(RichPlot3)

# Create an empty data frame that will hold the output data
RichPlotSumm = data.frame(Comp=levels(RichPlot3$Comp))
RichPlotSumm$Estimate = 0
RichPlotSumm$Error = 0
RichPlotSumm$p = 0
head(RichPlotSumm)

# Run a for loop that goes through each subset of data and makes the estimates using the betta function
for (i in levels(RichPlot3$Comp)){
  d = RichPlot3[RichPlot3$Comp==i,]
  Betta = betta(d$Estimate,d$Error)
  print(Betta$table)
  RichPlotSumm[RichPlotSumm$Comp==i,]$Estimate = Betta$table[,1]
  RichPlotSumm[RichPlotSumm$Comp==i,]$Error = Betta$table[,2]
  RichPlotSumm[RichPlotSumm$Comp==i,]$p = Betta$table[,3]
}

# Pull out sample type and facilty from the Comp variable so we can use them as plotting variables
lst <- strsplit(RichPlotSumm$Comp,' ')
v1 <- lapply(lst, `[`, 1)
v2 <- lapply(lst, `[`, 2)
RichPlotSumm$facility = v1
RichPlotSumm$sample.type = v2
RichPlotSumm

# Plot final richness estimates with error bars
# ±1.96*SE represents 95% confidence intervals
# Arlington plot first
RichPlotSummArlington = RichPlotSumm[RichPlotSumm$facility=="Arlington",]
RichPlotSummArlington
RichPlotSummArlington$sample.type = factor(RichPlotSummArlington$sample.type, levels = c("Endo","Ecto","Manure"))
p = ggplot(RichPlotSummArlington,aes(y=Estimate,x=sample.type))
p = p + geom_point(size=3) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), width=0.2)
p = p + theme_bw()
p = p + ylab("Richness estimate")
p = p + ggtitle("Arlington")
p

# Plot final richness estimates with error bars
# ±1.96*SE represents 95% confidence intervals
# Now DCC
RichPlotSummDCC = RichPlotSumm[RichPlotSumm$facility=="DCC",]
RichPlotSummDCC
RichPlotSummDCC$sample.type = factor(RichPlotSummDCC$sample.type, levels = c("Endo","Ecto","Manure"))
p = ggplot(RichPlotSummDCC,aes(y=Estimate,x=sample.type))
p = p + geom_point(size=3) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), width=0.2)
p = p + theme_bw()
p = p + ylab("Richness estimate")
p = p + ggtitle("DCC")
p
