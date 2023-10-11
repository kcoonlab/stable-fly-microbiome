set.seed(123456)

source("./Helper_Functions.R") #From CITATION

# Load required packages
library(ANCOMBC)
library(breakaway)
library(DescTools)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(phyloseq)
library(purrr)
library(tidyverse)
library(vegan)

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

# Calculate median read count by sample type
sdt = data.frame(as(sample_data(ps), "data.frame"),
                 TotalReads = sample_sums(ps), keep.rownames = TRUE)
MySummary <- sdt %>%
  group_by(sample.type) %>%
  summarize(median_count = median(TotalReads, na.rm=TRUE)) 
MySummary %>% print(n = Inf)

# Generate breakaway richness estimates

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

#### Goal: Summarize breakaway estimates across different sample types (internal flies, external flies, manure) within each facility (Arlington, DCC) ####

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

# Plot final richness estimates with 95% confidence intervals

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

# Now DCC plot
RichPlotSummDCC = RichPlotSumm[RichPlotSumm$facility=="DCC",]
RichPlotSummDCC
RichPlotSummDCC$sample.type = factor(RichPlotSummDCC$sample.type, levels = c("Endo","Ecto","Manure"))
p = ggplot(RichPlotSummDCC,aes(y=Estimate,x=sample.type))
p = p + geom_point(size=3) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), width=0.2)
p = p + theme_bw()
p = p + ylab("Richness estimate")
p = p + ggtitle("DCC")
p

#### Goal: Summarize the estimates across different sample types/sampling dates (Arlington samples only) ####

# First, make a single variable that contains all of these elements
RichPlot3Arlington = RichPlot3[RichPlot3$location=="Arlington",]
RichPlot3Arlington$Comp = paste(RichPlot3Arlington$sampling.date,RichPlot3Arlington$sample.type)
RichPlot3Arlington$Comp = as.factor(RichPlot3Arlington$Comp)
head(RichPlot3Arlington)

# Create an empty data frame that will hold the output data
RichPlotSummArlington = data.frame(Comp=levels(RichPlot3Arlington$Comp))
RichPlotSummArlington$Estimate = 0
RichPlotSummArlington$Error = 0
RichPlotSummArlington$p = 0
head(RichPlotSummArlington)

# Run a for loop that goes through each subset of data and makes the estimates using the betta function
for (i in levels(RichPlot3Arlington$Comp)){
  d = RichPlot3Arlington[RichPlot3Arlington$Comp==i,]
  Betta = betta(d$Estimate,d$Error)
  print(Betta$table)
  RichPlotSummArlington[RichPlotSummArlington$Comp==i,]$Estimate = Betta$table[,1]
  RichPlotSummArlington[RichPlotSummArlington$Comp==i,]$Error = Betta$table[,2]
  RichPlotSummArlington[RichPlotSummArlington$Comp==i,]$p = Betta$table[,3]
}

# Pull out elements from the Comp variable so we can use them as plotting variables
lst <- strsplit(RichPlotSummArlington$Comp,' ')
v1 <- lapply(lst, `[`, 1)
v2 <- lapply(lst, `[`, 2)
RichPlotSummArlington$sampling.date = v1
RichPlotSummArlington$sample.type = v2
RichPlotSummArlington

# Plot final richness estimates with 95% confidence intervals
RichPlotSummArlington$sampling.date = factor(RichPlotSummArlington$sampling.date, levels = c("7.9.2021","7.16.2021","7.23.2021","7.30.2021","8.6.2021","8.13.2021","8.20.2021","8.27.2021","9.10.2021"))
RichPlotSummArlington$sample.type = factor(RichPlotSummArlington$sample.type, levels = c("Endo","Ecto","Manure"))
p = ggplot(RichPlotSummArlington,aes(y=Estimate,x=sampling.date,color=sample.type))
p = p + geom_point(size=3) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), width=0.2)
p = p + theme_bw()
p = p + ylab("Richness estimate")
p = p + ggtitle("Arlington")
p

#### Goal: Summarize the estimates across different sample types/sampling dates (DCC samples only) ####

# Start with manure samples
RichPlot3DCCManure = RichPlot3[RichPlot3$location=="DCC"&RichPlot3$sample.type=="Manure",]
RichPlot3DCCManure$sampling.date = as.factor(RichPlot3DCCManure$sampling.date)
head(RichPlot3DCCManure)

# Create an empty data frame that will hold the output data
RichPlotSummDCCManure = data.frame(sampling.date=levels(RichPlot3DCCManure$sampling.date))
RichPlotSummDCCManure$Estimate = 0
RichPlotSummDCCManure$Error = 0
RichPlotSummDCCManure$p = 0
head(RichPlotSummDCCManure)

# Run a for loop that goes through each subset of data and makes the estimates using the betta function
for (i in levels(RichPlot3DCCManure$sampling.date)){
  d = RichPlot3DCCManure[RichPlot3DCCManure$sampling.date==i,]
  Betta = betta(d$Estimate,d$Error)
  print(Betta$table)
  RichPlotSummDCCManure[RichPlotSummDCCManure$sampling.date==i,]$Estimate = Betta$table[,1]
  RichPlotSummDCCManure[RichPlotSummDCCManure$sampling.date==i,]$Error = Betta$table[,2]
  RichPlotSummDCCManure[RichPlotSummDCCManure$sampling.date==i,]$p = Betta$table[,3]
}

# Plot final richness estimates with 95% confidence intervals
RichPlotSummDCCManure$sampling.date = factor(RichPlotSummDCCManure$sampling.date, levels = c("7.8.2021","7.15.2021","7.22.2021","7.29.2021","8.5.2021","8.12.2021","8.19.2021","8.26.2021","9.2.2021","9.9.2021","9.15.2021"))
p = ggplot(RichPlotSummDCCManure,aes(y=Estimate,x=sampling.date))
p = p + geom_point(size=3) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), width=0.2)
p = p + theme_bw()
p = p + ylab("Richness estimate")
p = p + ggtitle("DCC")
p

# Now endo samples (for dates having at least three independent pools)
RichPlot3DCCEndo = RichPlot3[RichPlot3$location=="DCC"&RichPlot3$sample.type=="Endo",]
RichPlot3DCCEndo = RichPlot3DCCEndo[RichPlot3DCCEndo$sampling.date=="7.15.2021"|RichPlot3DCCEndo$sampling.date=="7.22.2021",]
RichPlot3DCCEndo$sampling.date = as.factor(RichPlot3DCCEndo$sampling.date)
head(RichPlot3DCCEndo)
RichPlotSummDCCEndo = data.frame(sampling.date=levels(RichPlot3DCCEndo$sampling.date))
RichPlotSummDCCEndo$Estimate = 0
RichPlotSummDCCEndo$Error = 0
RichPlotSummDCCEndo$p = 0
head(RichPlotSummDCCEndo)
for (i in levels(RichPlot3DCCEndo$sampling.date)){
  d = RichPlot3DCCEndo[RichPlot3DCCEndo$sampling.date==i,]
  Betta = betta(d$Estimate,d$Error)
  print(Betta$table)
  RichPlotSummDCCEndo[RichPlotSummDCCEndo$sampling.date==i,]$Estimate = Betta$table[,1]
  RichPlotSummDCCEndo[RichPlotSummDCCEndo$sampling.date==i,]$Error = Betta$table[,2]
  RichPlotSummDCCEndo[RichPlotSummDCCEndo$sampling.date==i,]$p = Betta$table[,3]
}
RichPlotSummDCCEndo$sampling.date = factor(RichPlotSummDCCEndo$sampling.date, levels = c("7.15.2021","7.22.2021"))
p = ggplot(RichPlotSummDCCEndo,aes(y=Estimate,x=sampling.date))
p = p + geom_point(size=3) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), width=0.2)
p = p + theme_bw()
p = p + ylab("Richness estimate")
p = p + ggtitle("DCC")
p

# Now endo samples (for dates having less than three independent pools)
RichPlot3DCCEndo = RichPlot3[RichPlot3$location=="DCC"&RichPlot3$sample.type=="Endo",]
RichPlot3DCCEndo = RichPlot3DCCEndo[RichPlot3DCCEndo$sampling.date!="7.15.2021",]
RichPlot3DCCEndo = RichPlot3DCCEndo[RichPlot3DCCEndo$sampling.date!="7.22.2021",]
RichPlot3DCCEndo$sampling.date = as.factor(RichPlot3DCCEndo$sampling.date)
RichPlot3DCCEndo
x=c(1:5) #7.8.2021,7.29.2021,8.5.2021,9.9.2021,9.15.2021
y=c(21.25279,33.07883,181.07717,73.04119,36.04553)
plot(x,y)

# Now ecto samples (all dates have less than three independent pools)
RichPlot3DCCEcto = RichPlot3[RichPlot3$location=="DCC"&RichPlot3$sample.type=="Ecto",]
RichPlot3DCCEcto$sampling.date = as.factor(RichPlot3DCCEcto$sampling.date)
RichPlot3DCCEcto
x=c(1:4) #7.22.2021,7.29.2021,8.5.2021,9.9.2021
y=c(19.09134,14.02042,66.40300,20.26297)
plot(x,y)

#### Goal: Summarize the estimates for each sample type across different sampling locations within a single facility ####

# Arlington manure first
RichPlot3ArlingtonManure = RichPlot3[RichPlot3$location=="Arlington"&RichPlot3$sample.type=="Manure",]
RichPlot3ArlingtonManure$trap = as.factor(RichPlot3ArlingtonManure$trap)
head(RichPlot3ArlingtonManure)

# Create an empty data frame that will hold the output data
RichPlotSummArlingtonManure = data.frame(trap=levels(RichPlot3ArlingtonManure$trap))
RichPlotSummArlingtonManure$Estimate = 0
RichPlotSummArlingtonManure$Error = 0
RichPlotSummArlingtonManure$p = 0
head(RichPlotSummArlingtonManure)

# Run a for loop that goes through each subset of data and makes the estimates using the betta function
for (i in levels(RichPlot3ArlingtonManure$trap)){
  d = RichPlot3ArlingtonManure[RichPlot3ArlingtonManure$trap==i,]
  Betta = betta(d$Estimate,d$Error)
  print(Betta$table)
  RichPlotSummArlingtonManure[RichPlotSummArlingtonManure$trap==i,]$Estimate = Betta$table[,1]
  RichPlotSummArlingtonManure[RichPlotSummArlingtonManure$trap==i,]$Error = Betta$table[,2]
  RichPlotSummArlingtonManure[RichPlotSummArlingtonManure$trap==i,]$p = Betta$table[,3]
}

# Plot final richness estimates with 95% confidence intervals
RichPlotSummArlingtonManure$trap = factor(RichPlotSummArlingtonManure$trap, levels = c("Arl-M1","Arl-M2","Arl-M3","Arl-M4","Arl-M5","Arl-sickpen"))
p = ggplot(RichPlotSummArlingtonManure,aes(y=Estimate,x=trap))
p = p + geom_point(size=3) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), width=0.2)
p = p + theme_bw()
p = p + ylab("Richness estimate")
p = p + ggtitle("Arlington Manure")
p

# Now Arlington flies (internal samples)
RichPlot3ArlingtonFlies = RichPlot3[RichPlot3$location=="Arlington"&RichPlot3$sample.type!="Manure",]
RichPlot3ArlingtonFlies = RichPlot3ArlingtonFlies[RichPlot3ArlingtonFlies$sample.type=="Endo",]
RichPlot3ArlingtonFlies$trap = as.factor(RichPlot3ArlingtonFlies$trap)
head(RichPlot3ArlingtonFlies)

# Create an empty data frame that will hold the output data
RichPlotSummArlingtonFlies = data.frame(trap=levels(RichPlot3ArlingtonFlies$trap))
RichPlotSummArlingtonFlies$Estimate = 0
RichPlotSummArlingtonFlies$Error = 0
RichPlotSummArlingtonFlies$p = 0
head(RichPlotSummArlingtonFlies)

# Run a for loop that goes through each subset of data and makes the estimates using the betta function
for (i in levels(RichPlot3ArlingtonFlies$trap)){
  d = RichPlot3ArlingtonFlies[RichPlot3ArlingtonFlies$trap==i,]
  Betta = betta(d$Estimate,d$Error)
  print(Betta$table)
  RichPlotSummArlingtonFlies[RichPlotSummArlingtonFlies$trap==i,]$Estimate = Betta$table[,1]
  RichPlotSummArlingtonFlies[RichPlotSummArlingtonFlies$trap==i,]$Error = Betta$table[,2]
  RichPlotSummArlingtonFlies[RichPlotSummArlingtonFlies$trap==i,]$p = Betta$table[,3]
}

# Plot final richness estimates with 95% confidence intervals
RichPlotSummArlingtonFlies$trap = factor(RichPlotSummArlingtonFlies$trap, levels = c("B1-East","B1-West","B2-East","B2-West"))
p = ggplot(RichPlotSummArlingtonFlies,aes(y=Estimate,x=trap))
p = p + geom_point(size=3) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), width=0.2)
p = p + theme_bw()
p = p + ylab("Richness estimate")
p = p + ggtitle("Arlington Flies (Internal)")
p

# Now Arlington flies (external samples)
RichPlot3ArlingtonFlies = RichPlot3[RichPlot3$location=="Arlington"&RichPlot3$sample.type!="Manure",]
RichPlot3ArlingtonFlies = RichPlot3ArlingtonFlies[RichPlot3ArlingtonFlies$sample.type=="Ecto",]
RichPlot3ArlingtonFlies$trap = as.factor(RichPlot3ArlingtonFlies$trap)
head(RichPlot3ArlingtonFlies)

# Create an empty data frame that will hold the output data
RichPlotSummArlingtonFlies = data.frame(trap=levels(RichPlot3ArlingtonFlies$trap))
RichPlotSummArlingtonFlies$Estimate = 0
RichPlotSummArlingtonFlies$Error = 0
RichPlotSummArlingtonFlies$p = 0
head(RichPlotSummArlingtonFlies)

# Run a for loop that goes through each subset of data and makes the estimates using the betta function
for (i in levels(RichPlot3ArlingtonFlies$trap)){
  d = RichPlot3ArlingtonFlies[RichPlot3ArlingtonFlies$trap==i,]
  Betta = betta(d$Estimate,d$Error)
  print(Betta$table)
  RichPlotSummArlingtonFlies[RichPlotSummArlingtonFlies$trap==i,]$Estimate = Betta$table[,1]
  RichPlotSummArlingtonFlies[RichPlotSummArlingtonFlies$trap==i,]$Error = Betta$table[,2]
  RichPlotSummArlingtonFlies[RichPlotSummArlingtonFlies$trap==i,]$p = Betta$table[,3]
}

# Plot final richness estimates with 95% confidence intervals
RichPlotSummArlingtonFlies$trap = factor(RichPlotSummArlingtonFlies$trap, levels = c("B1-East","B1-West","B2-East","B2-West"))
p = ggplot(RichPlotSummArlingtonFlies,aes(y=Estimate,x=trap))
p = p + geom_point(size=3) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), width=0.2)
p = p + theme_bw()
p = p + ylab("Richness estimate")
p = p + ggtitle("Arlington Flies (External)")
p

# Now DCC manure 
RichPlot3DCCManure = RichPlot3[RichPlot3$location=="DCC"&RichPlot3$sample.type=="Manure",]
RichPlot3DCCManure$trap = as.factor(RichPlot3DCCManure$trap)
head(RichPlot3DCCManure)

# Create an empty data frame that will hold the output data
RichPlotSummDCCManure = data.frame(trap=levels(RichPlot3DCCManure$trap))
RichPlotSummDCCManure$Estimate = 0
RichPlotSummDCCManure$Error = 0
RichPlotSummDCCManure$p = 0
head(RichPlotSummDCCManure)

# Run a for loop that goes through each subset of data and makes the estimates using the betta function
for (i in levels(RichPlot3DCCManure$trap)){
  d = RichPlot3DCCManure[RichPlot3DCCManure$trap==i,]
  Betta = betta(d$Estimate,d$Error)
  print(Betta$table)
  RichPlotSummDCCManure[RichPlotSummDCCManure$trap==i,]$Estimate = Betta$table[,1]
  RichPlotSummDCCManure[RichPlotSummDCCManure$trap==i,]$Error = Betta$table[,2]
  RichPlotSummDCCManure[RichPlotSummDCCManure$trap==i,]$p = Betta$table[,3]
}

# Plot final richness estimates with 95% confidence intervals
RichPlotSummDCCManure$trap = factor(RichPlotSummDCCManure$trap, levels = c("DCC-Q1","DCC-Q2","DCC-Q3","DCC-Q4","DCC-Outdoor"))
p = ggplot(RichPlotSummDCCManure,aes(y=Estimate,x=trap))
p = p + geom_point(size=3) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), width=0.2)
p = p + theme_bw()
p = p + ylab("Richness estimate")
p = p + ggtitle("DCC Manure")
p

# Calculate descriptive stats RE: bacterial diversity in different sample types

# Total ASVs in manure vs. internal fly vs. external fly samples
ps.Manure <- subset_samples(ps, sample.type == "Manure")
asvfac <- rownames(otu_table(ps.Manure))
asvtab = apply(otu_table(ps.Manure), MARGIN = 2, function(x) {
    tapply(x, INDEX = asvfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
sum(apply(asvtab > 1, 1, sum) != 0)  

ps.Endo <- subset_samples(ps, sample.type == "Endo")
asvfac <- rownames(otu_table(ps.Endo))
asvtab = apply(otu_table(ps.Endo), MARGIN = 2, function(x) {
    tapply(x, INDEX = asvfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
sum(apply(asvtab > 1, 1, sum) != 0)  

ps.Ecto <- subset_samples(ps, sample.type == "Ecto")
asvfac <- rownames(otu_table(ps.Ecto))
asvtab = apply(otu_table(ps.Ecto), MARGIN = 2, function(x) {
    tapply(x, INDEX = asvfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
sum(apply(asvtab > 1, 1, sum) != 0)  

# Total bacterial orders represented across all samples
ordfac = factor(tax_table(ps)[, "Order"])
ordtab = apply(otu_table(ps), MARGIN = 2, function(x) {
    tapply(x, INDEX = ordfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
nrow(ordtab)

#Calculate mean relative abundance of abundant bacterial orders (>2% in at least one sample type on a given sampling date) in different sample types from both facilities
ps.Order <- tax_glom(ps, taxrank = "Order", NArm = FALSE)
ps.Order.sub <- subset_taxa(ps.Order, Order %in% c(" Acetobacterales"," Caulobacterales"," Rhizobiales"," Rhodobacterales"," Sphingomonadales",
									" Aeromonadales"," Burkholderiales"," Cardiobacteriales"," Chromatiales"," Enterobacterales",
									" Pseudomonadales"," Xanthomonadales"," Bacillales"," Exiguobacterales"," Lactobacillales",
									" Staphylococcales"," Christensenellales"," Clostridiales"," Clostridia"," Eubacteriales",
									" Lachnospirales"," Monoglobales"," Oscillospirales"," Peptostreptococcales-Tissierellales"," Acholeplasmatales",
									" Erysipelotrichales"," Flavobacteriales"," Bacteroidales"," Bifidobacteriales"," Corynebacteriales",
									" Micrococcales"," Propionibacteriales"," Coriobacteriales"," Campylobacterales"," Spirochaetales"))
relabun.ps.Order.sub <- transform_sample_counts(ps.Order.sub,function(x) x / sum(x))
df.relabun.ps.Order.sub <- psmelt(relabun.ps.Order.sub)
MySummary <- df.relabun.ps.Order.sub %>%
  group_by(summary, Order) %>%
  summarize(mean_abund = mean(Abundance, na.rm=TRUE))
MySummary %>% print(n = Inf)

# Differential abundance analysis of microbiome counts between different sample types/locations

# Arlington vs. DCC-derived internal fly samples
ps.Order.sub.Endo = subset_samples(ps.Order.sub, summary=="Endo.Arlington" | summary=="Endo.DCC")
ancombc.out = ancombc(phyloseq = ps.Order.sub.Endo, formula = "summary", p_adj_method = "fdr")
df.ancombc.out <- data.frame(ancombc.out$res)
df.ancombc.out <- df.ancombc.out %>%
  filter(q_val.summaryEndo.DCC < 0.05) %>%
  arrange(lfc.summaryEndo.DCC)
df.ancombc.out

# Arlington vs. DCC-derived manure samples
ps.Order.sub.Manure = subset_samples(ps.Order.sub, pooltype=="Manure_Samples")
ancombc.out = ancombc(phyloseq = ps.Order.sub.Manure, formula = "location", p_adj_method = "fdr")
df.ancombc.out <- data.frame(ancombc.out$res) 
df.ancombc.out <- df.ancombc.out %>%
  filter(q_val.locationDCC < 0.05) %>%
  arrange(lfc.locationDCC)
df.ancombc.out

#### Taxonomic Bar Plots ####
#### Goal: Relative abundance of bacterial orders in different sample types (by sampling date/facility) ####

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
        
orderAbdDF.Formatted <- read.csv('orderAbdDF_formatted.csv',header=TRUE)
metadata <- metadata %>%
  rownames_to_column(var = "Sample")
sampleLevels <- as.character(metadata$Sample)
sampleLabels <- metadata %>%
  dplyr::select(Sample)
sampleLabelsVec <- sampleLabels$Sample
names(sampleLabelsVec) <- sampleLabels$Sample
orderAbdDF.Formatted$Sample <- factor(orderAbdDF.Formatted$Sample, levels = sampleLevels)
orderAbdDF.Formatted$Order <- factor(orderAbdDF.Formatted$Order, levels=c(" Acetobacterales"," Caulobacterales"," Rhizobiales"," Rhodobacterales"," Sphingomonadales",
									" Aeromonadales"," Burkholderiales"," Cardiobacteriales"," Chromatiales"," Enterobacterales",
									" Pseudomonadales"," Xanthomonadales"," Bacillales"," Exiguobacterales"," Lactobacillales",
									" Staphylococcales"," Christensenellales"," Clostridiales"," Clostridia"," Eubacteriales",
									" Lachnospirales"," Monoglobales"," Oscillospirales"," Peptostreptococcales-Tissierellales"," Acholeplasmatales",
									" Erysipelotrichales"," Flavobacteriales"," Bacteroidales"," Bifidobacteriales"," Corynebacteriales",
									" Micrococcales"," Propionibacteriales"," Coriobacteriales"," Campylobacterales"," Spirochaetales"))
cbPalette <- unique(orderAbdDF.Formatted$Color)
orderCompPlot <- ggplot(orderAbdDF.Formatted,
                         aes_string(x = "Sample", y = "Abundance",
                                    fill = "Order")) +
  geom_bar(stat = "identity", width = 1) +
  theme_bw() +
  scale_x_discrete(labels = sampleLabelsVec) +
  scale_fill_manual(values=cbPalette)
orderCompPlot

#### Goal: PCoA to compare communities in different sample types from each facility ####

# Arlington samples first
                                                
# Normalize ps to relative abundances
ps.norm = transform_sample_counts(ps, function(x) x / sum(x) )
                                  
# rename for ease of use
ps = ps.norm
                                  
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

# Now DCC samples
				  
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

#### Goal: PCoA to compare communities in different sample types across sampling dates within each facility ####

# Arlington samples first
ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="7.9.2021")
ps.ordination.PCoA = ordinate(ps.arlington, method="PCoA", distance="bray")
p1 = plot_ordination(ps.arlington, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p1 = p1 + theme_bw() + geom_point(size =1.5)
p1 = p1 + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34", #Endo
                                      "#19468C", #Ecto
                                      "#B75463"),#Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p1 = p1 + theme(legend.text = element_text(size = 10))
p1$layers = p1$layers[-1] #to remove the larger point coded in original plot
p1

ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="7.16.2021")
ps.ordination.PCoA = ordinate(ps.arlington, method="PCoA", distance="bray")
p2 = plot_ordination(ps.arlington, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p2 = p2 + theme_bw() + geom_point(size =1.5)
p2 = p2 + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34", #Endo
                                      "#19468C", #Ecto
                                      "#B75463"),#Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p2 = p2 + theme(legend.text = element_text(size = 10))
p2$layers = p2$layers[-1] #to remove the larger point coded in original plot
p2
				  
ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="7.23.2021")
ps.ordination.PCoA = ordinate(ps.arlington, method="PCoA", distance="bray")
p3 = plot_ordination(ps.arlington, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p3 = p3 + theme_bw() + geom_point(size =1.5)
p3 = p3 + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34", #Endo
                                      "#19468C", #Ecto
                                      "#B75463"),#Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p3 = p3 + theme(legend.text = element_text(size = 10))
p3$layers = p3$layers[-1] #to remove the larger point coded in original plot
p3
				  
ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="7.30.2021")
ps.ordination.PCoA = ordinate(ps.arlington, method="PCoA", distance="bray")
p4 = plot_ordination(ps.arlington, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p4 = p4 + theme_bw() + geom_point(size =1.5)
p4 = p4 + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34", #Endo
                                      "#19468C", #Ecto
                                      "#B75463"),#Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p4 = p4 + theme(legend.text = element_text(size = 10))
p4$layers = p4$layers[-1] #to remove the larger point coded in original plot
p4
				  
ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="8.6.2021")
ps.ordination.PCoA = ordinate(ps.arlington, method="PCoA", distance="bray")
p5 = plot_ordination(ps.arlington, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p5 = p5 + theme_bw() + geom_point(size =1.5)
p5 = p5 + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34", #Endo
                                      "#19468C", #Ecto
                                      "#B75463"),#Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p5 = p5 + theme(legend.text = element_text(size = 10))
p5$layers = p5$layers[-1] #to remove the larger point coded in original plot
p5
				  
ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="8.13.2021")
ps.ordination.PCoA = ordinate(ps.arlington, method="PCoA", distance="bray")
p6 = plot_ordination(ps.arlington, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p6 = p6 + theme_bw() + geom_point(size =1.5)
p6 = p6 + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34", #Endo
                                      "#19468C", #Ecto
                                      "#B75463"),#Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p6 = p6 + theme(legend.text = element_text(size = 10))
p6$layers = p6$layers[-1] #to remove the larger point coded in original plot
p6
				  
ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="8.20.2021")
ps.ordination.PCoA = ordinate(ps.arlington, method="PCoA", distance="bray")
p7 = plot_ordination(ps.arlington, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p7 = p7 + theme_bw() + geom_point(size =1.5)
p7 = p7 + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34", #Endo
                                      "#19468C", #Ecto
                                      "#B75463"),#Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p7 = p7 + theme(legend.text = element_text(size = 10))
p7$layers = p7$layers[-1] #to remove the larger point coded in original plot
p7
				  
ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="8.27.2021")
ps.ordination.PCoA = ordinate(ps.arlington, method="PCoA", distance="bray")
p8 = plot_ordination(ps.arlington, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p8 = p8 + theme_bw() + geom_point(size =1.5)
p8 = p8 + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34", #Endo
                                      "#19468C", #Ecto
                                      "#B75463"),#Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p8 = p8 + theme(legend.text = element_text(size = 10))
p8$layers = p8$layers[-1] #to remove the larger point coded in original plot
p8
				  
ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="9.10.2021")
ps.ordination.PCoA = ordinate(ps.arlington, method="PCoA", distance="bray")
p9 = plot_ordination(ps.arlington, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p9 = p9 + theme_bw() + geom_point(size =1.5)
p9 = p9 + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34",  #Endo
                                      "#19468C",  #Ecto
                                      "#B75463"), #Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p9 = p9 + theme(legend.text = element_text(size = 10))
p9$layers = p9$layers[-1] #to remove the larger point coded in original plot
p9

ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9 + rremove("x.text"), 
          labels = c("Jul9", "Jul16", "Jul23", "Jul30", "Aug6", "Aug13", "Aug20", "Aug27", "Sept10"),
          ncol = 3, nrow = 3)

# Now DCC samples
ps.dcc <- subset_samples(ps, location=="DCC"&sampling.date=="7.8.2021")
ps.ordination.PCoA = ordinate(ps.dcc, method="PCoA", distance="bray")
p1 = plot_ordination(ps.dcc, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p1 = p1 + theme_bw() + geom_point(size =1.5)
p1 = p1 + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34", #Endo
                                      "#19468C", #Ecto
                                      "#B75463"),#Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p1 = p1 + theme(legend.text = element_text(size = 10))
p1$layers = p1$layers[-1] #to remove the larger point coded in original plot
p1

ps.dcc <- subset_samples(ps, location=="DCC"&sampling.date=="7.15.2021")
ps.ordination.PCoA = ordinate(ps.dcc, method="PCoA", distance="bray")
p2 = plot_ordination(ps.dcc, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p2 = p2 + theme_bw() + geom_point(size =1.5)
p2 = p2 + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34", #Endo
                                      "#19468C", #Ecto
                                      "#B75463"),#Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p2 = p2 + theme(legend.text = element_text(size = 10))
p2$layers = p2$layers[-1] #to remove the larger point coded in original plot
p2
				  
ps.dcc <- subset_samples(ps, location=="DCC"&sampling.date=="7.22.2021")
ps.ordination.PCoA = ordinate(ps.dcc, method="PCoA", distance="bray")
p3 = plot_ordination(ps.dcc, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p3 = p3 + theme_bw() + geom_point(size =1.5)
p3 = p3 + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34", #Endo
                                      "#19468C", #Ecto
                                      "#B75463"),#Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p3 = p3 + theme(legend.text = element_text(size = 10))
p3$layers = p3$layers[-1] #to remove the larger point coded in original plot
p3
				  
ps.dcc <- subset_samples(ps, location=="DCC"&sampling.date=="7.29.2021")
ps.ordination.PCoA = ordinate(ps.dcc, method="PCoA", distance="bray")
p4 = plot_ordination(ps.dcc, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p4 = p4 + theme_bw() + geom_point(size =1.5)
p4 = p4 + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34", #Endo
                                      "#19468C", #Ecto
                                      "#B75463"),#Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p4 = p4 + theme(legend.text = element_text(size = 10))
p4$layers = p4$layers[-1] #to remove the larger point coded in original plot
p4
				  
ps.dcc <- subset_samples(ps, location=="DCC"&sampling.date=="8.5.2021")
ps.ordination.PCoA = ordinate(ps.dcc, method="PCoA", distance="bray")
p5 = plot_ordination(ps.dcc, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p5 = p5 + theme_bw() + geom_point(size =1.5)
p5 = p5 + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34", #Endo
                                      "#19468C", #Ecto
                                      "#B75463"),#Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p5 = p5 + theme(legend.text = element_text(size = 10))
p5$layers = p5$layers[-1] #to remove the larger point coded in original plot
p5
				  
ps.dcc <- subset_samples(ps, location=="DCC"&sampling.date=="9.9.2021")
ps.ordination.PCoA = ordinate(ps.dcc, method="PCoA", distance="bray")
p6 = plot_ordination(ps.dcc, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p6 = p6 + theme_bw() + geom_point(size =1.5)
p6 = p6 + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34", #Endo
                                      "#19468C", #Ecto
                                      "#B75463"),#Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p6 = p6 + theme(legend.text = element_text(size = 10))
p6$layers = p6$layers[-1] #to remove the larger point coded in original plot
p6
				  
ps.dcc <- subset_samples(ps, location=="DCC"&sampling.date=="9.15.2021")
ps.ordination.PCoA = ordinate(ps.dcc, method="PCoA", distance="bray")
p7 = plot_ordination(ps.dcc, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "sample.type", label = NULL, justDF = FALSE)
p7 = p7 + theme_bw() + geom_point(size =1.5)
p7 = p7 + scale_color_manual(breaks = c("Endo", "Ecto", "Manure"),
                           values = c("#DFAB34", #Endo
                                      "#19468C", #Ecto
                                      "#B75463"),#Manure
                           name="", 
                           labels=c("Internal","External","Manure"))
p7 = p7 + theme(legend.text = element_text(size = 10))
p7$layers = p7$layers[-1] #to remove the larger point coded in original plot
p7

ggarrange(p1, p2, p3, p4, p5, p6, p7 + rremove("x.text"), 
          labels = c("Jul8", "Jul15", "Jul22", "Jul29", "Aug5", "Sept9", "Sept15"),
          ncol = 3, nrow = 3)

#### Goal: PCoA to compare communities in different sample types across different sampling locations within a single facility ####

# Arlington manure first
                                  
# Create PCoA ordination
ps.arlington <- subset_samples(ps, location=="Arlington"&sample.type=="Manure")
ps.ordination.PCoA = ordinate(ps.arlington, method="PCoA", distance="bray")

# plot ordination
p = plot_ordination(ps.arlington, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "trap", label = NULL, justDF = FALSE)
p = p + theme_bw() + geom_point(size =1.5)
p = p + theme(legend.text = element_text(size = 10))
p$layers = p$layers[-1] #to remove the larger point coded in original plot
p

# Now Arlington flies

# Create PCoA ordination
ps.arlington <- subset_samples(ps, location=="Arlington"&sample.type!="Manure")
ps.ordination.PCoA = ordinate(ps.arlington, method="PCoA", distance="bray")

# plot ordination
p = plot_ordination(ps.arlington, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "trap", label = NULL, justDF = FALSE)
p = p + theme_bw() + geom_point(size =1.5)
p = p + theme(legend.text = element_text(size = 10))
p$layers = p$layers[-1] #to remove the larger point coded in original plot
p

# Now DCC manure

# Create PCoA ordination
ps.dcc <- subset_samples(ps, location=="DCC"&sample.type=="Manure")
ps.ordination.PCoA = ordinate(ps.dcc, method="PCoA", distance="bray")

# plot ordination
p = plot_ordination(ps.dcc, ps.ordination.PCoA, type = "samples", 
                axes = 1:2, color = "trap", label = NULL, justDF = FALSE)
p = p + theme_bw() + geom_point(size =1.5)
p = p + theme(legend.text = element_text(size = 10))
p$layers = p$layers[-1] #to remove the larger point coded in original plot
p

#### PERMANOVA on BC-Dissimilarity ####
#### Goal: We want to know whether the distances between samples correspond to their source ####

# Arlington samples first

# Create a variable that is your distance matrix; You'll need this later.
ps.arlington <- subset_samples(ps, location=="Arlington")
DistVar = phyloseq::distance(ps.arlington, method = "bray")

# Extract the sample_data from the phyloseq object and turn it into a dataframe; You'll need this later.
psData = data.frame(sample_data(ps.arlington))

# Run PERMANOVA by sample type
adonis2(DistVar ~ sample.type, data = psData, method = "bray") # P = 0.001, R2 = 0.20856

# Create a function for pairwise comparisons, to test significant effects between each treatment
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='BH',reduce=NULL,perm=999)
{
  co <- combn(unique(as.character(factors)),2)
  pairs <- c()
  Df <- c()
  SumsOfSqs <- c()
  F.Model <- c()
  R2 <- c()
  p.value <- c()
  
  for(elem in 1:ncol(co)){
    if(inherits(x, 'dist')){
      x1=as.matrix(x)[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                      factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
    }
    
    else  (
      if (sim.function == 'daisy'){
        x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
      } 
      else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    )
    
    ad <- adonis2(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])],
                 permutations = perm);
    pairs <- c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    Df <- c(Df,ad[1,1])
    SumsOfSqs <- c(SumsOfSqs, ad[1,2])
    F.Model <- c(F.Model,ad[1,4]);
    R2 <- c(R2,ad[1,3]);
    p.value <- c(p.value,ad[1,5])
  }
  p.adjusted <- p.adjust(p.value,method=p.adjust.m)
  
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  pairw.res <- data.frame(pairs,Df,SumsOfSqs,F.Model,R2,p.value,p.adjusted,sig)
  
  if(!is.null(reduce)){
    pairw.res <- subset (pairw.res, grepl(reduce,pairs))
    pairw.res$p.adjusted <- p.adjust(pairw.res$p.value,method=p.adjust.m)
    
    sig = c(rep('',length(pairw.res$p.adjusted)))
    sig[pairw.res$p.adjusted <= 0.1] <-'.'
    sig[pairw.res$p.adjusted <= 0.05] <-'*'
    sig[pairw.res$p.adjusted <= 0.01] <-'**'
    sig[pairw.res$p.adjusted <= 0.001] <-'***'
    pairw.res <- data.frame(pairw.res[,1:7],sig)
  }
  class(pairw.res) <- c("pwadonis", "data.frame")
  return(pairw.res)
} 

# Run pairwise adonis function on distance variable, if there was a significant overall effect
pairwise.adonis(DistVar,psData$sample.type) # Endo vs Manure, P = 0.001, R2 = 0.21925948
					    # Endo vs Ecto, P = 0.001, R2 = 0.02240995
					    # Manure vs Ecto, P = 0.001, R2 = 0.27261527

# Statistical comparison of disperson within each sample type
beta = betadisper(DistVar, psData$sample.type)
anova(beta) # P < 0.0001
test = permutest(beta, pairwise=TRUE, permutations=999)
test
Tukey=TukeyHSD(beta, conf.level = 0.95)
Tukey   # Endo-Ecto, p = 0.4666097
	# Manure-Ecto, p < 0.0001
	# Manure-Endo, p < 0.0001

# Now DCC samples
ps.dcc <- subset_samples(ps, location=="DCC")
DistVar = phyloseq::distance(ps.dcc, method = "bray")
psData = data.frame(sample_data(ps.dcc))
adonis2(DistVar ~ sample.type, data = psData, method = "bray") # P = 0.001, R2 = 0.27593
				  
pairwise.adonis(DistVar,psData$sample.type) # Endo vs Manure, P = 0.0015, R2 = 0.23744811
					    # Endo vs Ecto, P = 0.2560, R2 = 0.06742574
					    # Manure vs Ecto, P = 0.0015, R2 = 0.16048601
				  
beta = betadisper(DistVar, psData$sample.type)
anova(beta) # P < 0.0001
test = permutest(beta, pairwise=TRUE, permutations=999)
test
				  
Tukey=TukeyHSD(beta, conf.level = 0.95)
Tukey   # Endo-Ecto, P = 0.0805425
	# Manure-Ecto, P = 0.0612392
	# Manure-Endo, P < 0.0001

#### PERMANOVA on BC-Dissimilarity ####
#### Goal: We want to know whether the distances between samples correspond to their source, irrespective of sampling date/facility ####

# Arlington samples first
ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="7.9.2021")
DistVar = phyloseq::distance(ps.arlington, method = "bray")
psData = data.frame(sample_data(ps.arlington))
adonis2(DistVar ~ sample.type, data = psData, method = "bray") # P = 0.001, R2 = 0.2914

beta = betadisper(DistVar, psData$sample.type)
anova(beta) # P = 0.008171
test = permutest(beta, pairwise=TRUE, permutations=999)
test # P = 0.004

ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="7.16.2021")
DistVar = phyloseq::distance(ps.arlington, method = "bray")
psData = data.frame(sample_data(ps.arlington))
adonis2(DistVar ~ sample.type, data = psData, method = "bray") # P = 0.001, R2 = 0.24975

co <- combn(unique(as.character(psData$sample.type)),2)  
x1=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,1]),as.character(co[2,1])),
	psData$sample.type %in% c(as.character(co[1,1]),as.character(co[2,1]))]    
x2=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,2]),as.character(co[2,2])),
	psData$sample.type %in% c(as.character(co[1,2]),as.character(co[2,2]))]    
x3=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,3]),as.character(co[2,3])),
	psData$sample.type %in% c(as.character(co[1,3]),as.character(co[2,3]))]    

ad1 <- adonis2(x1 ~ psData$sample.type[psData$sample.type %in% c(co[1,1],co[2,1])], permutations = 999)
ad2 <- adonis2(x2 ~ psData$sample.type[psData$sample.type %in% c(co[1,2],co[2,2])], permutations = 999)
ad3 <- adonis2(x3 ~ psData$sample.type[psData$sample.type %in% c(co[1,3],co[2,3])], permutations = 999)

R2 <- c(ad1[1,3],ad2[1,3],ad3[1,3])
p.values <- c(ad1[1,5],ad2[1,5],ad3[1,5])
p.adjusted <- p.adjust(p.values,method='BH')
R2
p.adjusted	# Endo-Manure, P = 0.003, R2 = 0.24749841
		# Endo-Ecto, P = 0.942
		# Manure-Ecto, P = 0.030, R2 = 0.41576891
  
beta = betadisper(DistVar, psData$sample.type)
anova(beta) # P = 0.0002489
test = permutest(beta, pairwise=TRUE, permutations=999)
test # P = 0.002
				  
Tukey=TukeyHSD(beta, conf.level = 0.95)
Tukey   # Endo-Ecto, P = 0.7426287
	# Manure-Ecto, P = 0.0189643
	# Manure-Endo, P = 0.0001726

ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="7.23.2021")
DistVar = phyloseq::distance(ps.arlington, method = "bray")
psData = data.frame(sample_data(ps.arlington))
adonis2(DistVar ~ sample.type, data = psData, method = "bray") # P = 0.001, R2 = 0.23475

co <- combn(unique(as.character(psData$sample.type)),2)  
x1=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,1]),as.character(co[2,1])),
	psData$sample.type %in% c(as.character(co[1,1]),as.character(co[2,1]))]    
x2=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,2]),as.character(co[2,2])),
	psData$sample.type %in% c(as.character(co[1,2]),as.character(co[2,2]))]    
x3=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,3]),as.character(co[2,3])),
	psData$sample.type %in% c(as.character(co[1,3]),as.character(co[2,3]))]    

ad1 <- adonis2(x1 ~ psData$sample.type[psData$sample.type %in% c(co[1,1],co[2,1])], permutations = 999)
ad2 <- adonis2(x2 ~ psData$sample.type[psData$sample.type %in% c(co[1,2],co[2,2])], permutations = 999)
ad3 <- adonis2(x3 ~ psData$sample.type[psData$sample.type %in% c(co[1,3],co[2,3])], permutations = 999)

R2 <- c(ad1[1,3],ad2[1,3],ad3[1,3])
p.values <- c(ad1[1,5],ad2[1,5],ad3[1,5])
p.adjusted <- p.adjust(p.values,method='BH') 
R2
p.adjusted	# Endo-Manure, P = 0.003, R2 = 0.24749841
		# Endo-Ecto, P = 0.942
		# Manure-Ecto, P = 0.030, R2 = 0.41576891
				  
beta = betadisper(DistVar, psData$sample.type)
anova(beta) # P < 0.0001
test = permutest(beta, pairwise=TRUE, permutations=999)
test # P = 0.001
				  
Tukey=TukeyHSD(beta, conf.level = 0.95)
Tukey   # Endo-Ecto, P = 0.1044805
	# Manure-Ecto, P = 0.0000160
	# Manure-Endo, P < 0.0001

ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="7.30.2021")
DistVar = phyloseq::distance(ps.arlington, method = "bray")
psData = data.frame(sample_data(ps.arlington))
adonis2(DistVar ~ sample.type, data = psData, method = "bray") # P = 0.001, R2 = 0.31227

co <- combn(unique(as.character(psData$sample.type)),2)  
x1=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,1]),as.character(co[2,1])),
	psData$sample.type %in% c(as.character(co[1,1]),as.character(co[2,1]))]    
x2=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,2]),as.character(co[2,2])),
	psData$sample.type %in% c(as.character(co[1,2]),as.character(co[2,2]))]    
x3=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,3]),as.character(co[2,3])),
	psData$sample.type %in% c(as.character(co[1,3]),as.character(co[2,3]))]    

ad1 <- adonis2(x1 ~ psData$sample.type[psData$sample.type %in% c(co[1,1],co[2,1])], permutations = 999)
ad2 <- adonis2(x2 ~ psData$sample.type[psData$sample.type %in% c(co[1,2],co[2,2])], permutations = 999)
ad3 <- adonis2(x3 ~ psData$sample.type[psData$sample.type %in% c(co[1,3],co[2,3])], permutations = 999)

R2 <- c(ad1[1,3],ad2[1,3],ad3[1,3])
p.values <- c(ad1[1,5],ad2[1,5],ad3[1,5])
p.adjusted <- p.adjust(p.values,method='BH') 
R2
p.adjusted	# Endo-Manure, P = 0.006, R2 = 0.30885951
		# Endo-Ecto, P = 0.684
		# Manure-Ecto, P = 0.012, R2 = 0.36516052
				  
beta = betadisper(DistVar, psData$sample.type)
anova(beta) # P = 0.009521
test = permutest(beta, pairwise=TRUE, permutations=999)
test # P = 0.003
				  
Tukey=TukeyHSD(beta, conf.level = 0.95)
Tukey   # Endo-Ecto, P = 0.6798374
	# Manure-Ecto, P = 0.0961909
	# Manure-Endo, P = 0.0080492
				  
ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="8.6.2021")
DistVar = phyloseq::distance(ps.arlington, method = "bray")
psData = data.frame(sample_data(ps.arlington))
adonis2(DistVar ~ sample.type, data = psData, method = "bray") # P = 0.001, R2 = 0.27656
				  
co <- combn(unique(as.character(psData$sample.type)),2)  
x1=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,1]),as.character(co[2,1])),
	psData$sample.type %in% c(as.character(co[1,1]),as.character(co[2,1]))]    
x2=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,2]),as.character(co[2,2])),
	psData$sample.type %in% c(as.character(co[1,2]),as.character(co[2,2]))]    
x3=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,3]),as.character(co[2,3])),
	psData$sample.type %in% c(as.character(co[1,3]),as.character(co[2,3]))]    
				  
ad1 <- adonis2(x1 ~ psData$sample.type[psData$sample.type %in% c(co[1,1],co[2,1])], permutations = 999)
ad2 <- adonis2(x2 ~ psData$sample.type[psData$sample.type %in% c(co[1,2],co[2,2])], permutations = 999)
ad3 <- adonis2(x3 ~ psData$sample.type[psData$sample.type %in% c(co[1,3],co[2,3])], permutations = 999)
				  
R2 <- c(ad1[1,3],ad2[1,3],ad3[1,3])
p.values <- c(ad1[1,5],ad2[1,5],ad3[1,5])
p.adjusted <- p.adjust(p.values,method='BH')
R2
p.adjusted	# Endo-Manure, P = 0.003, R2 = 0.24231427
		# Endo-Ecto, P = 0.007, R2 = 0.09929417
		# Manure-Ecto, P = 0.007, R2 = 0.43754122
				  
beta = betadisper(DistVar, psData$sample.type)
anova(beta) # P = 0.001843
test = permutest(beta, pairwise=TRUE, permutations=999)
test # P = 0.004
				  
Tukey=TukeyHSD(beta, conf.level = 0.95)
Tukey   # Endo-Ecto, P = 0.1086115
	# Manure-Ecto, P = 0.6108099
	# Manure-Endo, P = 0.0022579
				  
ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="8.13.2021")
DistVar = phyloseq::distance(ps.arlington, method = "bray")
psData = data.frame(sample_data(ps.arlington))
adonis2(DistVar ~ sample.type, data = psData, method = "bray") # P = 0.001, R2 = 0.30899
				  
co <- combn(unique(as.character(psData$sample.type)),2)  
x1=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,1]),as.character(co[2,1])),
	psData$sample.type %in% c(as.character(co[1,1]),as.character(co[2,1]))]    
x2=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,2]),as.character(co[2,2])),
	psData$sample.type %in% c(as.character(co[1,2]),as.character(co[2,2]))]    
x3=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,3]),as.character(co[2,3])),
	psData$sample.type %in% c(as.character(co[1,3]),as.character(co[2,3]))]  
				  
ad1 <- adonis2(x1 ~ psData$sample.type[psData$sample.type %in% c(co[1,1],co[2,1])], permutations = 999)
ad2 <- adonis2(x2 ~ psData$sample.type[psData$sample.type %in% c(co[1,2],co[2,2])], permutations = 999)
ad3 <- adonis2(x3 ~ psData$sample.type[psData$sample.type %in% c(co[1,3],co[2,3])], permutations = 999)
				  
R2 <- c(ad1[1,3],ad2[1,3],ad3[1,3])
p.values <- c(ad1[1,5],ad2[1,5],ad3[1,5])
p.adjusted <- p.adjust(p.values,method='BH')
R2 
p.adjusted	# Endo-Manure, P = 0.003, R2 = 0.29691997
		# Endo-Ecto, P = 0.546
		# Manure-Ecto, P = 0.009, R2 = 0.44078437
				  
beta = betadisper(DistVar, psData$sample.type)
anova(beta) # P = 0.07576
test = permutest(beta, pairwise=TRUE, permutations=999)
test # P = 0.082
				  
ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="8.20.2021")
DistVar = phyloseq::distance(ps.arlington, method = "bray")
psData = data.frame(sample_data(ps.arlington))
adonis2(DistVar ~ sample.type, data = psData, method = "bray") # P = 0.001, R2 = 0.33157
				  
co <- combn(unique(as.character(psData$sample.type)),2)  
x1=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,1]),as.character(co[2,1])),
	psData$sample.type %in% c(as.character(co[1,1]),as.character(co[2,1]))]    
x2=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,2]),as.character(co[2,2])),
	psData$sample.type %in% c(as.character(co[1,2]),as.character(co[2,2]))]    
x3=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,3]),as.character(co[2,3])),
	psData$sample.type %in% c(as.character(co[1,3]),as.character(co[2,3]))]    
				  
ad1 <- adonis2(x1 ~ psData$sample.type[psData$sample.type %in% c(co[1,1],co[2,1])], permutations = 999)
ad2 <- adonis2(x2 ~ psData$sample.type[psData$sample.type %in% c(co[1,2],co[2,2])], permutations = 999)
ad3 <- adonis2(x3 ~ psData$sample.type[psData$sample.type %in% c(co[1,3],co[2,3])], permutations = 999)
				  
R2 <- c(ad1[1,3],ad2[1,3],ad3[1,3])
p.values <- c(ad1[1,5],ad2[1,5],ad3[1,5])
p.adjusted <- p.adjust(p.values,method='BH')
R2
p.adjusted	# Endo-Manure, P = 0.006, R2 = 0.3166126
		# Endo-Ecto, P = 0.091
		# Manure-Ecto, P = 0.009, R2 = 0.3643591
				  
beta = betadisper(DistVar, psData$sample.type)
anova(beta) # P = 0.001678
test = permutest(beta, pairwise=TRUE, permutations=999)
test # P = 0.001
				  
Tukey=TukeyHSD(beta, conf.level = 0.95)
Tukey   # Endo-Ecto, P = 0.3872561
	# Manure-Ecto, P = 0.0649469
	# Manure-Endo, P = 0.0012142

ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="8.27.2021")
DistVar = phyloseq::distance(ps.arlington, method = "bray")
psData = data.frame(sample_data(ps.arlington))
adonis2(DistVar ~ sample.type, data = psData, method = "bray") # P = 0.001, R2 = 0.38069
				  
co <- combn(unique(as.character(psData$sample.type)),2)  
x1=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,1]),as.character(co[2,1])),
	psData$sample.type %in% c(as.character(co[1,1]),as.character(co[2,1]))]    
x2=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,2]),as.character(co[2,2])),
	psData$sample.type %in% c(as.character(co[1,2]),as.character(co[2,2]))]    
x3=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,3]),as.character(co[2,3])),
	psData$sample.type %in% c(as.character(co[1,3]),as.character(co[2,3]))]    
				  
ad1 <- adonis2(x1 ~ psData$sample.type[psData$sample.type %in% c(co[1,1],co[2,1])], permutations = 999)
ad2 <- adonis2(x2 ~ psData$sample.type[psData$sample.type %in% c(co[1,2],co[2,2])], permutations = 999)
ad3 <- adonis2(x3 ~ psData$sample.type[psData$sample.type %in% c(co[1,3],co[2,3])], permutations = 999)
				  
R2 <- c(ad1[1,3],ad2[1,3],ad3[1,3])
p.values <- c(ad1[1,5],ad2[1,5],ad3[1,5])
p.adjusted <- p.adjust(p.values,method='BH')
R2
p.adjusted	# Endo-Manure, P = 0.003, R2 = 0.3456423
		# Endo-Ecto, P = 0.060
		# Manure-Ecto, P = 0.003, R2 = 0.5577478
				  
beta = betadisper(DistVar, psData$sample.type)
anova(beta) # P = 0.1226
test = permutest(beta, pairwise=TRUE, permutations=999)
test # P = 0.12

ps.arlington <- subset_samples(ps, location=="Arlington"&sampling.date=="9.10.2021")
DistVar = phyloseq::distance(ps.arlington, method = "bray")
psData = data.frame(sample_data(ps.arlington))
adonis2(DistVar ~ sample.type, data = psData, method = "bray") # P = 0.001, R2 = 0.33471
				  
co <- combn(unique(as.character(psData$sample.type)),2)  
x1=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,1]),as.character(co[2,1])),
	psData$sample.type %in% c(as.character(co[1,1]),as.character(co[2,1]))]    
x2=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,2]),as.character(co[2,2])),
	psData$sample.type %in% c(as.character(co[1,2]),as.character(co[2,2]))]    
x3=as.matrix(DistVar)[psData$sample.type %in% c(as.character(co[1,3]),as.character(co[2,3])),
	psData$sample.type %in% c(as.character(co[1,3]),as.character(co[2,3]))]    
				  
ad1 <- adonis2(x1 ~ psData$sample.type[psData$sample.type %in% c(co[1,1],co[2,1])], permutations = 999)
ad2 <- adonis2(x2 ~ psData$sample.type[psData$sample.type %in% c(co[1,2],co[2,2])], permutations = 999)
ad3 <- adonis2(x3 ~ psData$sample.type[psData$sample.type %in% c(co[1,3],co[2,3])], permutations = 999)
				  
R2 <- c(ad1[1,3],ad2[1,3],ad3[1,3])
p.values <- c(ad1[1,5],ad2[1,5],ad3[1,5])
p.adjusted <- p.adjust(p.values,method='BH')
R2
p.adjusted	# Endo-Manure, P = 0.003, R2 = 0.31185141
		# Endo-Ecto, P = 0.4020
		# Manure-Ecto, P = 0.0555, R2 = 0.37089031
				  
beta = betadisper(DistVar, psData$sample.type)
anova(beta) # P = 0.06986
test = permutest(beta, pairwise=TRUE, permutations=999)
test  P = 0.082

# Now DCC samples
ps.dcc <- subset_samples(ps, location=="DCC"&sampling.date=="7.15.2021")
DistVar = phyloseq::distance(ps.dcc, method = "bray")
psData = data.frame(sample_data(ps.dcc))
adonis2(DistVar ~ sample.type, data = psData, method = "bray") # P = 0.005, R2 = 0.42856
				  
beta = betadisper(DistVar, psData$sample.type)
anova(beta) # P = 0.000901
test = permutest(beta, pairwise=TRUE, permutations=999)
test # P = 0.001

ps.dcc <- subset_samples(ps, location=="DCC"&sampling.date=="7.22.2021")
DistVar = phyloseq::distance(ps.dcc, method = "bray")
psData = data.frame(sample_data(ps.dcc))
adonis2(DistVar ~ sample.type, data = psData, method = "bray") # P = 0.005, R2 = 0.40706
				  
beta = betadisper(DistVar, psData$sample.type)
anova(beta) # P = 0.0002146
test = permutest(beta, pairwise=TRUE, permutations=999)
test # P = 0.001

#### PERMANOVA on BC-Dissimilarity ####
#### We want to know whether the distances between samples correspond to their sampling location within a given facility ####

# Arlington manure samples first
ps.arlington <- subset_samples(ps, location=="Arlington"&sample.type=="Manure")
DistVar = phyloseq::distance(ps.arlington, method = "bray")
psData = data.frame(sample_data(ps.arlington))
adonis2(DistVar ~ trap, data = psData, method = "bray") # P = 0.001, R2 = 0.20768
				  
pairwise.adonis(DistVar,psData$trap) # Arl-M2 vs Arl-sickpen, P = 0.0075000
				     # Arl-M4 vs Arl-sickpen, P = 0.0100000
				     # Arl-M5 vs Arl-sickpen, P = 0.0075000
				  
beta = betadisper(DistVar, psData$trap)
anova(beta) # P = 0.1888
test = permutest(beta, pairwise=TRUE, permutations=999)
test # P = 0.168

# Now Arlington flies
ps.arlington <- subset_samples(ps, location=="Arlington"&sample.type!="Manure")
DistVar = phyloseq::distance(ps.arlington, method = "bray")
psData = data.frame(sample_data(ps.arlington))
adonis2(DistVar ~ trap, data = psData, method = "bray") # P = 0.007, R2 = 0.03325
pairwise.adonis(DistVar,psData$trap) # All NS (P > 0.05)

beta = betadisper(DistVar, psData$trap)
anova(beta) # P = 0.831
test = permutest(beta, pairwise=TRUE, permutations=999)
test # P = 0.833

# Now DCC manure
ps.dcc <- subset_samples(ps, location=="DCC"&sample.type=="Manure")
DistVar = phyloseq::distance(ps.dcc, method = "bray")
psData = data.frame(sample_data(ps.dcc))
adonis2(DistVar ~ trap, data = psData, method = "bray") # P = 0.004, R2 = 0.10697
pairwise.adonis(DistVar,psData$trap) # All NS (P > 0.05)
				  
beta = betadisper(DistVar, psData$trap)
anova(beta) # P = 0.2989
test = permutest(beta, pairwise=TRUE, permutations=999)
test # P = 0.285

#### Taxonomic Bar Plots ####
#### Goal: Relative abundance of bacterial orders in different sample type (by sampling location/facility) ####

ps.Maure <- subset_samples(ps, sample.type == "Manure")
variable1 = as.character(get_variable(ps.Maure, "sample.type"))
variable2 = as.character(get_variable(ps.Maure, "location"))
variable3 = as.character(get_variable(ps.Maure, "trap"))
sample_data(ps.Maure)$NewPastedVar <- mapply(paste0, variable1, variable2, variable3, collapse = "_")
ps.Merged <- merge_samples(ps.Maure, "NewPastedVar")
metadata <- data.frame(sample_data(ps.Merged))
orderAbdDF.Manure.byLocation <- MakeAbundanceDF(physeq = ps.Merged,
                               taxRank = "Order",
                               abundanceFilter = 0.02)
# write.csv(orderAbdDF.Manure.byLocation, "orderAbdDF.byLocation.csv", row.names=FALSE)
        
orderAbdDF.Manure.byLocation.Formatted <- read.csv('orderAbdDF.Manure.byLocation_formatted.csv',header=TRUE)
metadata <- metadata %>%
  rownames_to_column(var = "Sample")
sampleLevels <- as.character(metadata$Sample)
sampleLabels <- metadata %>%
  dplyr::select(Sample)
sampleLabelsVec <- sampleLabels$Sample
names(sampleLabelsVec) <- sampleLabels$Sample
orderAbdDF.Manure.byLocation.Formatted$Sample <- factor(orderAbdDF.Manure.byLocation.Formatted$Sample, levels = sampleLevels)
orderAbdDF.Manure.byLocation.Formatted$Order <- factor(orderAbdDF.Manure.byLocation.Formatted$Order, levels=c(" Burkholderiales"," Pseudomonadales"," Bacillales"," Lactobacillales",
						" Staphylococcales"," Christensenellales"," Lachnospirales"," Monoglobales"," Oscillospirales",
						" Peptostreptococcales-Tissierellales"," Acholeplasmatales"," Erysipelotrichales"," Flavobacteriales"," Bacteroidales",
						" Bifidobacteriales"," Corynebacteriales"," Spirochaetales"))
cbPalette <- unique(orderAbdDF.Manure.byLocation.Formatted$Color)
orderCompPlot <- ggplot(orderAbdDF.Manure.byLocation.Formatted,
                         aes_string(x = "Sample", y = "Abundance",
                                    fill = "Order")) +
  geom_bar(stat = "identity", width = 1) +
  theme_bw() +
  scale_x_discrete(labels = sampleLabelsVec) +
  scale_fill_manual(values=cbPalette)
orderCompPlot				  
