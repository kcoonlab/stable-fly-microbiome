set.seed(123456)

library(phyloseq)
library(dplyr)
library(tidyverse)
library(data.table)
library(breakaway)
library(ggpubr)

#### Fig. S1

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

#### Goal: Summarize the estimates across different sample types/sampling dates
#### Fig. S1 (left panel, Arlington samples only)

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

# Plot final richness estimates with error bars
# ±1.96*SE represents 95% confidence intervals
RichPlotSummArlington$sampling.date = factor(RichPlotSummArlington$sampling.date, levels = c("7.9.2021","7.16.2021","7.23.2021","7.30.2021","8.6.2021","8.13.2021","8.20.2021","8.27.2021","9.10.2021"))
RichPlotSummArlington$sample.type = factor(RichPlotSummArlington$sample.type, levels = c("Endo","Ecto","Manure"))
p = ggplot(RichPlotSummArlington,aes(y=Estimate,x=sampling.date,color=sample.type))
p = p + geom_point(size=3) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), width=0.2)
p = p + theme_bw()
p = p + ylab("Richness estimate")
p = p + ggtitle("Arlington")
p

#### Goal: Summarize the estimates across different sample types/sampling dates
#### Fig. S1 (right panel, DCC samples only)
#### Start with manure samples

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

# Plot final richness estimates with error bars
# ±1.96*SE represents 95% confidence intervals
RichPlotSummDCCManure$sampling.date = factor(RichPlotSummDCCManure$sampling.date, levels = c("7.8.2021","7.15.2021","7.22.2021","7.29.2021","8.5.2021","8.12.2021","8.19.2021","8.26.2021","9.2.2021","9.9.2021","9.15.2021"))
p = ggplot(RichPlotSummDCCManure,aes(y=Estimate,x=sampling.date))
p = p + geom_point(size=3) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), width=0.2)
p = p + theme_bw()
p = p + ylab("Richness estimate")
p = p + ggtitle("DCC")
p

#### Now endo samples (for dates having at least three independent pools)

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

#### Now endo samples (for dates having less than three independent pools)

RichPlot3DCCEndo = RichPlot3[RichPlot3$location=="DCC"&RichPlot3$sample.type=="Endo",]
RichPlot3DCCEndo = RichPlot3DCCEndo[RichPlot3DCCEndo$sampling.date!="7.15.2021",]
RichPlot3DCCEndo = RichPlot3DCCEndo[RichPlot3DCCEndo$sampling.date!="7.22.2021",]
RichPlot3DCCEndo$sampling.date = as.factor(RichPlot3DCCEndo$sampling.date)
RichPlot3DCCEndo
x=c(1:5) #7.8.2021,7.29.2021,8.5.2021,9.9.2021,9.15.2021
y=c(21.25279,33.07883,181.07717,73.04119,36.04553)
plot(x,y)

#### Now ecto samples (all dates have less than three independent pools)

RichPlot3DCCEcto = RichPlot3[RichPlot3$location=="DCC"&RichPlot3$sample.type=="Ecto",]
RichPlot3DCCEcto$sampling.date = as.factor(RichPlot3DCCEcto$sampling.date)
RichPlot3DCCEcto
x=c(1:4) #7.22.2021,7.29.2021,8.5.2021,9.9.2021
y=c(19.09134,14.02042,66.40300,20.26297)
plot(x,y)

#### Goal: summarize the estimates for each sample type across different sampling locations within a single facility
#### Fig. S2 (left panel, Arlington samples only)

#### Arlington manure first

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

# Plot final richness estimates with error bars
# ±1.96*SE represents 95% confidence intervals
RichPlotSummArlingtonManure$trap = factor(RichPlotSummArlingtonManure$trap, levels = c("Arl-M1","Arl-M2","Arl-M3","Arl-M4","Arl-M5","Arl-sickpen"))
p = ggplot(RichPlotSummArlingtonManure,aes(y=Estimate,x=trap))
p = p + geom_point(size=3) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), width=0.2)
p = p + theme_bw()
p = p + ylab("Richness estimate")
p = p + ggtitle("Arlington Manure")
p

#### Now Arlington flies (internal samples)

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

# Plot final richness estimates with error bars
# ±1.96*SE represents 95% confidence intervals
RichPlotSummArlingtonFlies$trap = factor(RichPlotSummArlingtonFlies$trap, levels = c("B1-East","B1-West","B2-East","B2-West"))
p = ggplot(RichPlotSummArlingtonFlies,aes(y=Estimate,x=trap))
p = p + geom_point(size=3) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), width=0.2)
p = p + theme_bw()
p = p + ylab("Richness estimate")
p = p + ggtitle("Arlington Flies (Internal)")
p

#### Now Arlington flies (external samples)

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

# Plot final richness estimates with error bars
# ±1.96*SE represents 95% confidence intervals
RichPlotSummArlingtonFlies$trap = factor(RichPlotSummArlingtonFlies$trap, levels = c("B1-East","B1-West","B2-East","B2-West"))
p = ggplot(RichPlotSummArlingtonFlies,aes(y=Estimate,x=trap))
p = p + geom_point(size=3) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), width=0.2)
p = p + theme_bw()
p = p + ylab("Richness estimate")
p = p + ggtitle("Arlington Flies (External)")
p

#### Goal: summarize the estimates for each sample type across different sampling locations within a single facility
#### Fig. S2 (right panel, DCC manure samples only)

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

# Plot final richness estimates with error bars
# ±1.96*SE represents 95% confidence intervals
RichPlotSummDCCManure$trap = factor(RichPlotSummDCCManure$trap, levels = c("DCC-Q1","DCC-Q2","DCC-Q3","DCC-Q4","DCC-Outdoor"))
p = ggplot(RichPlotSummDCCManure,aes(y=Estimate,x=trap))
p = p + geom_point(size=3) + geom_errorbar(aes(ymin=Estimate-1.96*Error,ymax=Estimate+1.96*Error), width=0.2)
p = p + theme_bw()
p = p + ylab("Richness estimate")
p = p + ggtitle("DCC Manure")
p

#### Goal: PCoA to compare communities in different sample types across sampling dates within each facility ####
#### Fig. S3 (Arlington samples only)

ps.norm = transform_sample_counts(ps, function(x) x / sum(x) )
# rename for ease of use.
ps = ps.norm

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

#### Goal: PCoA to compare communities in different sample types across sampling dates within each facility ####
#### Fig. S4 (DCC samples only)

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
#### Fig. S5 (top left panel, Arlington manure samples only)

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

#### Fig. S5 (bottom left panel, Arlington fly samples only)

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

#### Fig. S5 (top right panel, DCC fly samples only)

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

#### Figs. S7-S9 (Arlington panels)

ps <- readRDS("ps_FieldWork2021_AJS_Final.rds")
sample_data(ps)$merge_factor <- c("Endo10A-7.16.2021","Endo10B-7.16.2021","Endo10D-7.16.2021","Endo11A-7.22.2021","Endo11B-7.22.2021","Endo11D-7.22.2021","Endo11E-7.22.2021","Endo12B-7.23.2021","Endo12C-7.23.2021","Endo12D-7.23.2021","Endo12E-7.23.2021","Endo13A-7.23.2021","Endo13B-7.23.2021","Endo13C-7.23.2021","Endo13D-7.23.2021","Endo13E-7.23.2021","Endo14A-7.23.2021","Endo14B-7.23.2021","Endo14C-7.23.2021","Endo15A-7.23.2021","Endo15B-7.23.2021","Endo15C-7.23.2021","Endo15D-7.23.2021","Endo16A-7.29.2021","Endo17C-7.30.2021","Endo18A-7.30.2021","Endo18B-7.30.2021","Endo19A-7.30.2021","Endo19B-7.30.2021","Endo1A-7.8.2021","Endo1B-7.8.2021","Endo20A-7.30.2021","Endo20B-7.30.2021","Endo21A-8.5.2021","Endo22A-8.6.2021","Endo22B-8.6.2021","Endo22C-8.6.2021","Endo22D-8.6.2021","Endo22E-8.6.2021","Endo23A-8.6.2021","Endo23B-8.6.2021","Endo23C-8.6.2021","Endo23D-8.6.2021","Endo23E-8.6.2021","Endo24A-8.6.2021","Endo24B-8.6.2021","Endo24C-8.6.2021","Endo24D-8.6.2021","Endo24E-8.6.2021","Endo25A-8.6.2021","Endo25B-8.6.2021","Endo25C-8.6.2021","Endo25D-8.6.2021","Endo27A-8.13.2021","Endo27B-8.13.2021","Endo27C-8.13.2021","Endo27D-8.13.2021","Endo27E-8.13.2021","Endo28A-8.13.2021","Endo28B-8.13.2021","Endo29A-8.13.2021","Endo29B-8.13.2021","Endo2A-7.9.2021","Endo2B-7.9.2021","Endo30A-8.13.2021","Endo30B-8.13.2021","Endo30C-8.13.2021","Endo30D-8.13.2021","Endo32A-8.20.2021","Endo32D-8.20.2021","Endo33A-8.20.2021","Endo33B-8.20.2021","Endo34A-8.20.2021","Endo34B-8.20.2021","Endo34C-8.20.2021","Endo35A-8.20.2021","Endo37A-8.27.2021","Endo37B-8.27.2021","Endo37C-8.27.2021","Endo37D-8.27.2021","Endo37E-8.27.2021","Endo38A-8.27.2021","Endo39A-8.27.2021","Endo39B-8.27.2021","Endo39C-8.27.2021","Endo39D-8.27.2021","Endo3A-7.9.2021","Endo3C-7.9.2021","Endo3D-7.9.2021","Endo3E-7.9.2021","Endo40A-8.27.2021","Endo40B-8.27.2021","Endo40C-8.27.2021","Endo41A-9.9.2021","Endo42A-9.10.2021","Endo42B-9.10.2021","Endo42C-9.10.2021","Endo42D-9.10.2021","Endo42E-9.10.2021","Endo44A-9.10.2021","Endo44B-9.10.2021","Endo44C-9.10.2021","Endo44D-9.10.2021","Endo44E-9.10.2021","Endo46A-9.15.2021","Endo4A-7.9.2021","Endo4B-7.9.2021","Endo4C-7.9.2021","Endo4D-7.9.2021","Endo5A-7.9.2021","Endo5B-7.9.2021","Endo6A-7.15.2021","Endo6B-7.15.2021","Endo6C-7.15.2021","Endo6E-7.15.2021","Endo7B-7.16.2021","Endo7C-7.16.2021","Endo7D-7.16.2021","Endo7E-7.16.2021","Endo8A-7.16.2021","Endo8C-7.16.2021","Endo8D-7.16.2021","Endo8E-7.16.2021","Endo9A-7.16.2021","Endo9B-7.16.2021","Endo9C-7.16.2021","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Ecto24-8.6.2021","Ecto10-7.16.2021","Ecto11-7.22.2021","Ecto12-7.23.2021","Ecto13-7.23.2021","Ecto14-7.23.2021","Ecto16-7.29.2021","Ecto17-7.30.2021","Ecto18-7.30.2021","Ecto19-7.30.2021","Ecto20-7.30.2021","Ecto21-8.5.2021","Ecto22-8.6.2021","Ecto23-8.6.2021","Ecto25-8.6.2021","Ecto27-8.13.2021","Ecto28-8.13.2021","Ecto29-8.13.2021","Ecto30-8.13.2021","Ecto32-8.20.2021","Ecto33-8.20.2021","Ecto34-8.20.2021","Ecto35-8.20.2021","Ecto37-8.27.2021","Ecto38-8.27.2021","Ecto39-8.27.2021","Ecto40-8.27.2021","Ecto41-9.9.2021","Ecto42-9.10.2021","Ecto44-9.10.2021","Ecto7-7.16.2021","Ecto9-7.16.2021")
ps.arlington <- subset_samples(ps, location=="Arlington")
ps.Merged <- merge_samples(ps.arlington,"merge_factor")
row_names <- row.names(sample_data(ps.Merged))
sample_data(ps.Merged)$sample_id <- row_names
lst <- strsplit(row_names,'-')
v1 <- lapply(lst, `[`, 1)
v2 <- lapply(lst, `[`, 2)
sample_data(ps.Merged)$pooltype <- v1
sample_data(ps.Merged)$sampling.date <- v2

setClass("Family",
         representation = representation(family_group = "character",
                                         member_list = "list",
                                         family_sum = "numeric"),
         prototype = prototype(family_group = NA_character_,
                               member_list = list(),
                               family_sum = NA_real_))

setClass("FamilyMember", 
         representation = representation(family_group = "character",
                                         name = "character",
                                         sample_designation = "character",
                                         sample_sum = "numeric",
                                         ASV_counts = "data.frame"),
         prototype = prototype(family_group = NA_character_,
                               name = NA_character_,
                               sample_designation = NA_character_,
                               sample_sum = NA_real_))

CountASVs <- function(sampleName, abundanceTable) {
  sampleAbundance <- abundanceTable %>%
    dplyr::select(sampleName) %>%
    rownames_to_column(var = "ASV")
  sampleAbundanceDF <- filter(sampleAbundance, sampleAbundance[, sampleName] > 0)
  return(sampleAbundanceDF)
}

GetSampleSum <- function(sampleName, physeqObject) {
  named_sum <- sample_sums(physeqObject)[sampleName]
  numeric_sum <- unname(named_sum)
}

# Pull out sample data
sampleData <- data.frame(sample_data(ps.Merged))
sampleData$sampling.date <- unlist(sampleData$sampling.date)
sampleData$pooltype <- unlist(sampleData$pooltype)

# Pull out ASV count table
abdTable <- as.matrix(otu_table(ps.Merged))
if (!taxa_are_rows(abdTable)) {
  abdTable <- t(abdTable)
}
abdTable <- as.data.frame(abdTable)

# Populate Classes

# This list will hold the Family objects once completed
familyClassList <- list()

# Loop over the rows of the sample data data frame, building the objects by
#   familiy group number.  But we only want to build it once, so when we run
#   across the family group number again we don't want to re-build it.

for (i in 1:140) {
  # Get the current family group from the sample data table
  currentFamilyGroup <- sampleData[i, "sampling.date"]
  
  # Get all the current family names in familyClassList to check if 
  #   currentFamilyGroup has already been built:
  familyNames <- names(familyClassList)
  
  # If the current family group has not been built, build it:
  if (!(currentFamilyGroup %in% familyNames)){
    # Populate the family_group slot:
    # Create a list of sample names for all the samples in the current family group
    #   by subsetting the sample data df and extracting the sample names as a char
    #   vector.  Then, add to the memberList list that will be looped over to
    #   build each FamilyMember object for the current Family.
    members <- sampleData %>%
      filter(sampling.date == currentFamilyGroup) %>%
      dplyr::select(all_of("sample_id"))
    memberList <- as.list(members$sample_id)
    
    # Create an empty list to hold FamilyMember objects which will be added to the 
    #   Family object in the end.
    currentFamilyMembersList <- list()
    
    # Loop over the sample names in memberList and build a new FamilyMember object for
    #   each sample:
    for (sample in memberList) {
      newMember <- new("FamilyMember",
                       family_group = currentFamilyGroup,
                       name = sample, 
                       sample_designation = sampleData$pooltype[sampleData$sample_id == sample],
                       sample_sum = GetSampleSum(sample, ps.Merged),
                       ASV_counts = CountASVs(sample, abdTable))
    # Add this new FamilyMember object to currentFamilyMemberList:
    currentFamilyMembersList[[sample]] <- newMember
    }
    
    # Now that the member list is complete, get a total read sum to be put into the 
    #  Family object (sum of each family member's total read count).
    totalSum <- 0
    for (member in currentFamilyMembersList) {
      sampleSum <- member@sample_sum
      totalSum <- sum(totalSum, sampleSum)
    }
    
    # We now have all the components needed to create new Family object.
    newFamily <- new("Family", 
                     family_group = currentFamilyGroup,
                     member_list = currentFamilyMembersList,
                     family_sum = totalSum)
    # Add to list
    familyClassList[[currentFamilyGroup]] <- newFamily
  }
}

# Select "Valid" Families

#Families including a manure sample and at least one fly sample.

#----- Determine which family groups are valid -----# 

# A valid Family must have at least one manure sample and one fly sample:
# So a Family object should have a family_member list equal to length 3 or
# equal to length 2 where one of the sample_designation values is "Manure"

# Create a list to hold the valid Family objects:
validFamilies <- familyClassList

# Infant Perspective 

#----- Populate Table for Plotting Infants -----#

# This is for the "infant perspective" plot.

# Goal: generate a data frame holding all infants and all their ASVs (as rows). For each
#   ASV in the given infant, count the number of reads and
#   the percentage of the total reads in the infant that the ASV accounts for
#   ("infant_percent" column.)

allF1ASVTables <- data.frame()

for (k in 1:length(validFamilies)) {
  family <- validFamilies[[k]]
  familyMembers <- family@member_list

  # Put each FamilyMember object in its own variable
  currentFly <- NULL

  for (i in 1:length(familyMembers)) {
      currentFly <- c(currentFly,familyMembers[[i]])
    	}

  # Put infants in a list to loop over:
  currentF1List <- list(currentFly)
  
  # Remove any NULL elements in the list (in case there is only 1 infant)
  currentF1List <- plyr::compact(currentF1List)

  # Loop over the infant list and populate the "infant perspective" data frame:
  #   allInfantASVTables - for each infant samples, holds the ASV, 
  #     its percentage of the total reads in the infant, and whether it's shared with mother
  #   (Later) to get the total percent shared and not shared for each infant (for plotting)
  #     allInfantASVTables %>% group_by(sample, shared) %>% summarise(percent = sum(infants))
  
  for (i in 1:length(currentF1List[[1]])) {
  
    current_F1 <- currentF1List[[1]][[i]]
    #currentInfantName <- current_infant@name
    
    # build infant's identifier (FamGroup.1 or FamGroup.2)
    F1Designation <- current_F1@sample_designation
    F1Identifier <- paste(family@family_group, F1Designation, sep = "/")
	
	Manure_ASV_counts = CountASVs("Manure-Manure", abdTable)

    # merge current infant's ASV counts with mother's
    F1Table <- merge(current_F1@ASV_counts,
                         Manure_ASV_counts,
                         by = "ASV", all = TRUE)
    
    # Replace column names with the samples' sample_designations for easier
    #   plotting downstream:
    
      currentColName <- names(F1Table)[2]
      setnames(F1Table, old = currentColName, new = familyMembers[[currentColName]]@sample_designation)
    
    # Collapse table to figure out what is shared from infant's perspective -
    #   remove any ASVs that have a count of NA in the current infant.
    F1Table <- F1Table[!(is.na(F1Table[[F1Designation]])), ]
    
    # Add a column to indicate if each ASV (row) in the collapsed table is 
    #   shared with mom or not. If the ASV has a value of NA in mother column
    #   then it cannot be shared between current infant and mom.
    
    F1Table$shared <- ifelse(!(is.na(F1Table$"Manure-Manure")), yes = TRUE, no = FALSE)
    F1Table$sample <- F1Identifier

    # Prepare to add this table to a bigger list, but don't add until end!
    setnames(F1Table, old = F1Designation, new = "read_count")
      
    # For this infant, calculate percent of reads that are from ASVs that
    #   are shared with mom, and not shared with mom (keep this table).
    F1Sum <- sum(F1Table$read_count)
    F1SharedSum <- sum(F1Table$read_count[F1Table$shared == TRUE])
    F1NotSharedSum <- sum(F1Table$read_count[F1Table$shared == FALSE])
      
    # Add the infantTable (the one with ASVs as rows) to larger data frame
    F1TablePercentages <- F1Table %>%
      mutate("F1_percent" = read_count/F1Sum) %>%
      dplyr::select(-c("Manure-Manure"))
    
    # Add infant sample name to table
    #infantTablePercentages <- cbind(infantTablePercentages,
                                    #"infant_sample" = currentInfantName)
    
    allF1ASVTables <- rbind(allF1ASVTables, F1TablePercentages)
  }
}

# Infant Perspective Richness

# Generate a plot that is similar to the one above, but by richness instead
#   of relative abundance of reads.

# Group allInfantASVTables by sample and shared (T/F) then count

F1PerspectiveRichness <- allF1ASVTables
array <- unlist(strsplit(F1PerspectiveRichness$sample, "[/]"))
date <- array[c(TRUE, FALSE)]
F1PerspectiveRichness$date <- date

F1PerspectiveRichness.1 <- F1PerspectiveRichness %>%
  group_by(date, sample, shared) %>%
  summarise(number_of_ASVs = n())
  
lst <- strsplit(F1PerspectiveRichness.1$sample,'/')
v2 <- lapply(lst, `[`, 2)
v2 = substr(v2,1,4)
F1PerspectiveRichness.1$sample.type <- v2

F1SharedOTUs <- F1PerspectiveRichness.1
F1SharedOTUs <- filter(F1SharedOTUs, shared == TRUE)
F1SharedOTUs <- subset(F1SharedOTUs, select = c(date,sample,number_of_ASVs, sample.type))
F1SharedOTUs$sample.type <- factor(F1SharedOTUs$sample.type, levels=c("Endo", "Ecto"))
F1SharedOTUs$trap <- c("B2-West","B1-East","B2-East","B2-West","B2-West","B2-West","B1-East","B1-East","B1-East","B1-East","B1-West","B1-West","B1-West","B1-West","B2-East","B2-East","B2-East","B1-East","B1-West","B2-East","B1-East","B1-East","B1-East","B1-East","B1-West","B1-West","B1-West","B1-West","B1-West","B2-East","B2-East","B2-East","B2-West","B2-West","B2-West","B2-West","B1-East","B1-West","B2-East","B2-West","B1-East","B1-West","B1-West","B2-East","B2-East","B2-West","B2-West","B1-East","B1-East","B1-West","B1-West","B1-West","B1-West","B2-East","B2-East","B2-East","B2-East","B2-West","B2-West","B1-East","B1-West","B2-East","B2-West","B1-East","B1-East","B1-East","B1-East","B1-East","B1-West","B1-West","B2-East","B2-East","B2-West","B2-West","B2-West","B2-West","B1-East","B1-West","B2-East","B2-West","B1-East","B1-East","B1-West","B1-West","B2-East","B2-East","B2-East","B2-West","B1-East","B1-West","B2-East","B2-West","B1-East","B1-East","B1-East","B1-East","B1-East","B1-West","B2-East","B2-East","B2-East","B2-East","B2-West","B2-West","B2-West","B1-East","B1-West","B2-East","B2-West","B1-East","B1-East","B1-East","B1-East","B1-East","B1-West","B1-West","B1-West","B1-West","B1-West","B2-East","B2-East","B2-East","B2-East","B2-East","B2-West","B2-West","B2-West","B2-West","B1-East","B2-East","B1-East","B1-East","B1-East","B1-East","B1-East","B2-East","B2-East","B2-East","B2-East","B2-East")
F1SharedOTUs$trap <- factor(F1SharedOTUs$trap, levels=c("B1-East", "B1-West", "B2-East", "B2-West"))
kruskal.test(number_of_ASVs~as.factor(sample.type),data=F1SharedOTUs) #Significant (p < 0.0001)
means <- aggregate(number_of_ASVs ~ sample.type, data = F1SharedOTUs, 
          FUN = function(x) c(mean = mean(x)))

cis <- aggregate(number_of_ASVs ~ sample.type, data = F1SharedOTUs, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))

summary <- merge(means,cis,by="sample.type")
summary$sample.type <- factor(summary$sample.type, levels=c("Endo", "Ecto"))

#### Fig. S7 (left, endo)                 
F1SharedOTUs$date <- factor(F1SharedOTUs$date, levels=c("7.9.2021","7.16.2021","7.23.2021","7.30.2021","8.6.2021","8.13.2021","8.20.2021","8.27.2021","9.10.2021"))
F1SharedOTUs.Endo <- F1SharedOTUs[F1SharedOTUs$sample.type=="Endo",]
F1SharedOTUs.Ecto <- F1SharedOTUs[F1SharedOTUs$sample.type=="Ecto",]
kruskal.test(number_of_ASVs~as.factor(date),data=F1SharedOTUs.Endo) #NS (p = 0.5275)
kruskal.test(number_of_ASVs~as.factor(date),data=F1SharedOTUs.Ecto) #NS (p = 0.0781)

means <- aggregate(number_of_ASVs ~ date, data = F1SharedOTUs.Endo, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(number_of_ASVs ~ date, data = F1SharedOTUs.Endo, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="date")
summary$date <- factor(summary$date, levels=c("7.9.2021","7.16.2021","7.23.2021","7.30.2021","8.6.2021","8.13.2021","8.20.2021","8.27.2021","9.10.2021"))
ggplot(summary, aes(x = date, y = number_of_ASVs.x)) +
    geom_errorbar(aes(ymin=number_of_ASVs.x-number_of_ASVs.y, ymax=number_of_ASVs.x+number_of_ASVs.y), width=.1) + #Fig. S7 (left)
    geom_point()

#### Fig. S7 (left, ecto)
means <- aggregate(number_of_ASVs ~ date, data = F1SharedOTUs.Ecto, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(number_of_ASVs ~ date, data = F1SharedOTUs.Ecto, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="date")
summary$date <- factor(summary$date, levels=c("7.9.2021","7.16.2021","7.23.2021","7.30.2021","8.6.2021","8.13.2021","8.20.2021","8.27.2021","9.10.2021"))
ggplot(summary, aes(x = date, y = number_of_ASVs.x)) +
    geom_errorbar(aes(ymin=number_of_ASVs.x-number_of_ASVs.y, ymax=number_of_ASVs.x+number_of_ASVs.y), width=.1) + #Fig. S7 (left)
    geom_point()

#### Fig. S8 (bottom, endo)
means <- aggregate(number_of_ASVs ~ trap, data = F1SharedOTUs.Endo, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(number_of_ASVs ~ trap, data = F1SharedOTUs.Endo, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="trap")
summary$date <- factor(summary$trap, levels=c("B1-East", "B1-West", "B2-East", "B2-West"))
ggplot(summary, aes(x = trap, y = number_of_ASVs.x)) + #Fig. S8 (bottom)
    geom_errorbar(aes(ymin=number_of_ASVs.x-number_of_ASVs.y, ymax=number_of_ASVs.x+number_of_ASVs.y), width=.1) +
    geom_point()

#### Fig. S8 (top, ecto)
means <- aggregate(number_of_ASVs ~ trap, data = F1SharedOTUs.Ecto, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(number_of_ASVs ~ trap, data = F1SharedOTUs.Ecto, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="trap")
summary$date <- factor(summary$trap, levels=c("B1-East", "B1-West", "B2-East", "B2-West"))
ggplot(summary, aes(x = date, y = number_of_ASVs.x)) + #Fig. S8 (top)
    geom_errorbar(aes(ymin=number_of_ASVs.x-number_of_ASVs.y, ymax=number_of_ASVs.x+number_of_ASVs.y), width=.1) +
    geom_point()

##

F1PerspectiveRA <- allF1ASVTables
array <- unlist(strsplit(F1PerspectiveRA$sample, "[/]"))
date <- array[c(TRUE, FALSE)]
F1PerspectiveRA$date <- date

F1PerspectiveRA.1 <- F1PerspectiveRA %>%
  group_by(date, sample, shared) %>%
  summarise(tot_shared_reads = sum(read_count))
  
F1PerspectiveRA.2 <- F1PerspectiveRA.1 %>%
  group_by(sample) %>%
  summarise(tot_reads = sum(tot_shared_reads))

F1PerspectiveRA.3 <- merge(x = F1PerspectiveRA.1,y = F1PerspectiveRA.2,by.x = "sample",by.y = "sample", all = T)

F1PerspectiveRA.3$percent <- F1PerspectiveRA.3$tot_shared_reads/F1PerspectiveRA.3$tot_reads

write.csv(F1PerspectiveRA.3, 
            file = "ManurePerspectiveRAFlies.csv",
            quote = FALSE, row.names = FALSE)

#### Fig. S9 (top left, endo)
data=read.csv("ManurePerspectiveRAFliesFormatted.csv")
data.endo <- data[data$fly.source=="Endo",]
data_summary <- aggregate(percent_shared ~ as.factor(date), data.endo,       # Create summary data
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("7.9.2021","7.16.2021","7.23.2021","7.30.2021","8.6.2021","8.13.2021","8.20.2021","8.27.2021","9.10.2021"))
ggplot(data_summary) +
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)

#### Fig. S9 (bottom left, ecto)
                          
data.ecto <- data[data$fly.source=="Ecto",]
data_summary <- aggregate(percent_shared ~ as.factor(date), data.ecto,       # Create summary data
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("7.9.2021","7.16.2021","7.23.2021","7.30.2021","8.6.2021","8.13.2021","8.20.2021","8.27.2021","9.10.2021"))
ggplot(data_summary) + #Fig. S9
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)









                          
#### Figs. S7-S9 (DCC panels)

ps <- readRDS("ps_FieldWork2021_AJS_Final.rds")
sample_data(ps)$merge_factor <- c("Endo10A-7.16.2021","Endo10B-7.16.2021","Endo10D-7.16.2021","Endo11A-7.22.2021","Endo11B-7.22.2021","Endo11D-7.22.2021","Endo11E-7.22.2021","Endo12B-7.23.2021","Endo12C-7.23.2021","Endo12D-7.23.2021","Endo12E-7.23.2021","Endo13A-7.23.2021","Endo13B-7.23.2021","Endo13C-7.23.2021","Endo13D-7.23.2021","Endo13E-7.23.2021","Endo14A-7.23.2021","Endo14B-7.23.2021","Endo14C-7.23.2021","Endo15A-7.23.2021","Endo15B-7.23.2021","Endo15C-7.23.2021","Endo15D-7.23.2021","Endo16A-7.29.2021","Endo17C-7.30.2021","Endo18A-7.30.2021","Endo18B-7.30.2021","Endo19A-7.30.2021","Endo19B-7.30.2021","Endo1A-7.8.2021","Endo1B-7.8.2021","Endo20A-7.30.2021","Endo20B-7.30.2021","Endo21A-8.5.2021","Endo22A-8.6.2021","Endo22B-8.6.2021","Endo22C-8.6.2021","Endo22D-8.6.2021","Endo22E-8.6.2021","Endo23A-8.6.2021","Endo23B-8.6.2021","Endo23C-8.6.2021","Endo23D-8.6.2021","Endo23E-8.6.2021","Endo24A-8.6.2021","Endo24B-8.6.2021","Endo24C-8.6.2021","Endo24D-8.6.2021","Endo24E-8.6.2021","Endo25A-8.6.2021","Endo25B-8.6.2021","Endo25C-8.6.2021","Endo25D-8.6.2021","Endo27A-8.13.2021","Endo27B-8.13.2021","Endo27C-8.13.2021","Endo27D-8.13.2021","Endo27E-8.13.2021","Endo28A-8.13.2021","Endo28B-8.13.2021","Endo29A-8.13.2021","Endo29B-8.13.2021","Endo2A-7.9.2021","Endo2B-7.9.2021","Endo30A-8.13.2021","Endo30B-8.13.2021","Endo30C-8.13.2021","Endo30D-8.13.2021","Endo32A-8.20.2021","Endo32D-8.20.2021","Endo33A-8.20.2021","Endo33B-8.20.2021","Endo34A-8.20.2021","Endo34B-8.20.2021","Endo34C-8.20.2021","Endo35A-8.20.2021","Endo37A-8.27.2021","Endo37B-8.27.2021","Endo37C-8.27.2021","Endo37D-8.27.2021","Endo37E-8.27.2021","Endo38A-8.27.2021","Endo39A-8.27.2021","Endo39B-8.27.2021","Endo39C-8.27.2021","Endo39D-8.27.2021","Endo3A-7.9.2021","Endo3C-7.9.2021","Endo3D-7.9.2021","Endo3E-7.9.2021","Endo40A-8.27.2021","Endo40B-8.27.2021","Endo40C-8.27.2021","Endo41A-9.9.2021","Endo42A-9.10.2021","Endo42B-9.10.2021","Endo42C-9.10.2021","Endo42D-9.10.2021","Endo42E-9.10.2021","Endo44A-9.10.2021","Endo44B-9.10.2021","Endo44C-9.10.2021","Endo44D-9.10.2021","Endo44E-9.10.2021","Endo46A-9.15.2021","Endo4A-7.9.2021","Endo4B-7.9.2021","Endo4C-7.9.2021","Endo4D-7.9.2021","Endo5A-7.9.2021","Endo5B-7.9.2021","Endo6A-7.15.2021","Endo6B-7.15.2021","Endo6C-7.15.2021","Endo6E-7.15.2021","Endo7B-7.16.2021","Endo7C-7.16.2021","Endo7D-7.16.2021","Endo7E-7.16.2021","Endo8A-7.16.2021","Endo8C-7.16.2021","Endo8D-7.16.2021","Endo8E-7.16.2021","Endo9A-7.16.2021","Endo9B-7.16.2021","Endo9C-7.16.2021","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Ecto24-8.6.2021","Ecto10-7.16.2021","Ecto11-7.22.2021","Ecto12-7.23.2021","Ecto13-7.23.2021","Ecto14-7.23.2021","Ecto16-7.29.2021","Ecto17-7.30.2021","Ecto18-7.30.2021","Ecto19-7.30.2021","Ecto20-7.30.2021","Ecto21-8.5.2021","Ecto22-8.6.2021","Ecto23-8.6.2021","Ecto25-8.6.2021","Ecto27-8.13.2021","Ecto28-8.13.2021","Ecto29-8.13.2021","Ecto30-8.13.2021","Ecto32-8.20.2021","Ecto33-8.20.2021","Ecto34-8.20.2021","Ecto35-8.20.2021","Ecto37-8.27.2021","Ecto38-8.27.2021","Ecto39-8.27.2021","Ecto40-8.27.2021","Ecto41-9.9.2021","Ecto42-9.10.2021","Ecto44-9.10.2021","Ecto7-7.16.2021","Ecto9-7.16.2021")
ps.dcc <- subset_samples(ps, location=="DCC")
ps.Merged <- merge_samples(ps.dcc,"merge_factor")
row_names <- row.names(sample_data(ps.Merged))
sample_data(ps.Merged)$sample_id <- row_names
lst <- strsplit(row_names,'-')
v1 <- lapply(lst, `[`, 1)
v2 <- lapply(lst, `[`, 2)
sample_data(ps.Merged)$pooltype <- v1
sample_data(ps.Merged)$sampling.date <- v2

setClass("Family",
         representation = representation(family_group = "character",
                                         member_list = "list",
                                         family_sum = "numeric"),
         prototype = prototype(family_group = NA_character_,
                               member_list = list(),
                               family_sum = NA_real_))

setClass("FamilyMember", 
         representation = representation(family_group = "character",
                                         name = "character",
                                         sample_designation = "character",
                                         sample_sum = "numeric",
                                         ASV_counts = "data.frame"),
         prototype = prototype(family_group = NA_character_,
                               name = NA_character_,
                               sample_designation = NA_character_,
                               sample_sum = NA_real_))

CountASVs <- function(sampleName, abundanceTable) {
  sampleAbundance <- abundanceTable %>%
    dplyr::select(sampleName) %>%
    rownames_to_column(var = "ASV")
  sampleAbundanceDF <- filter(sampleAbundance, sampleAbundance[, sampleName] > 0)
  return(sampleAbundanceDF)
}

GetSampleSum <- function(sampleName, physeqObject) {
  named_sum <- sample_sums(physeqObject)[sampleName]
  numeric_sum <- unname(named_sum)
}

# Pull out sample data
sampleData <- data.frame(sample_data(ps.Merged))
sampleData$sampling.date <- unlist(sampleData$sampling.date)
sampleData$pooltype <- unlist(sampleData$pooltype)

# Pull out ASV count table
abdTable <- as.matrix(otu_table(ps.Merged))
if (!taxa_are_rows(abdTable)) {
  abdTable <- t(abdTable)
}
abdTable <- as.data.frame(abdTable)

# Populate Classes

# This list will hold the Family objects once completed
familyClassList <- list()

# Loop over the rows of the sample data data frame, building the objects by
#   familiy group number.  But we only want to build it once, so when we run
#   across the family group number again we don't want to re-build it.

for (i in 1:18) {
  # Get the current family group from the sample data table
  currentFamilyGroup <- sampleData[i, "sampling.date"]
  
  # Get all the current family names in familyClassList to check if 
  #   currentFamilyGroup has already been built:
  familyNames <- names(familyClassList)
  
  # If the current family group has not been built, build it:
  if (!(currentFamilyGroup %in% familyNames)){
    # Populate the family_group slot:
    # Create a list of sample names for all the samples in the current family group
    #   by subsetting the sample data df and extracting the sample names as a char
    #   vector.  Then, add to the memberList list that will be looped over to
    #   build each FamilyMember object for the current Family.
    members <- sampleData %>%
      filter(sampling.date == currentFamilyGroup) %>%
      dplyr::select(all_of("sample_id"))
    memberList <- as.list(members$sample_id)
    
    # Create an empty list to hold FamilyMember objects which will be added to the 
    #   Family object in the end.
    currentFamilyMembersList <- list()
    
    # Loop over the sample names in memberList and build a new FamilyMember object for
    #   each sample:
    for (sample in memberList) {
      newMember <- new("FamilyMember",
                       family_group = currentFamilyGroup,
                       name = sample, 
                       sample_designation = sampleData$pooltype[sampleData$sample_id == sample],
                       sample_sum = GetSampleSum(sample, ps.Merged),
                       ASV_counts = CountASVs(sample, abdTable))
    # Add this new FamilyMember object to currentFamilyMemberList:
    currentFamilyMembersList[[sample]] <- newMember
    }
    
    # Now that the member list is complete, get a total read sum to be put into the 
    #  Family object (sum of each family member's total read count).
    totalSum <- 0
    for (member in currentFamilyMembersList) {
      sampleSum <- member@sample_sum
      totalSum <- sum(totalSum, sampleSum)
    }
    
    # We now have all the components needed to create new Family object.
    newFamily <- new("Family", 
                     family_group = currentFamilyGroup,
                     member_list = currentFamilyMembersList,
                     family_sum = totalSum)
    # Add to list
    familyClassList[[currentFamilyGroup]] <- newFamily
  }
}

# Select "Valid" Families

#Families including a mother and at least one infant.

#----- Determine which family groups are valid -----# 

# A valid Family must have at least one mother and one infant:
# So a Family object should have a family_member list equal to length 3 or
# equal to length 2 where one of the sample_designation values is "mother"

# Create a list to hold the valid Family objects:
validFamilies <- familyClassList

# Infant Perspective 

#----- Populate Table for Plotting Infants -----#

# This is for the "infant perspective" plot.

# Goal: generate a data frame holding all infants and all their ASVs (as rows). For each
#   ASV in the given infant, count the number of reads and
#   the percentage of the total reads in the infant that the ASV accounts for
#   ("infant_percent" column.)

allF1ASVTables <- data.frame()

for (k in 1:length(validFamilies)) {
  family <- validFamilies[[k]]
  familyMembers <- family@member_list

  # Put each FamilyMember object in its own variable
  currentFly <- NULL

  for (i in 1:length(familyMembers)) {
      currentFly <- c(currentFly,familyMembers[[i]])
    	}

  # Put infants in a list to loop over:
  currentF1List <- list(currentFly)
  
  # Remove any NULL elements in the list (in case there is only 1 infant)
  currentF1List <- plyr::compact(currentF1List)

  # Loop over the infant list and populate the "infant perspective" data frame:
  #   allInfantASVTables - for each infant samples, holds the ASV, 
  #     its percentage of the total reads in the infant, and whether it's shared with mother
  #   (Later) to get the total percent shared and not shared for each infant (for plotting)
  #     allInfantASVTables %>% group_by(sample, shared) %>% summarise(percent = sum(infants))
  
  for (i in 1:length(currentF1List[[1]])) {
  
    current_F1 <- currentF1List[[1]][[i]]
    #currentInfantName <- current_infant@name
    
    # build infant's identifier (FamGroup.1 or FamGroup.2)
    F1Designation <- current_F1@sample_designation
    F1Identifier <- paste(family@family_group, F1Designation, sep = "/")
	
	Manure_ASV_counts = CountASVs("Manure-Manure", abdTable)

    # merge current infant's ASV counts with mother's
    F1Table <- merge(current_F1@ASV_counts,
                         Manure_ASV_counts,
                         by = "ASV", all = TRUE)
    
    # Replace column names with the samples' sample_designations for easier
    #   plotting downstream:
    
      currentColName <- names(F1Table)[2]
      setnames(F1Table, old = currentColName, new = familyMembers[[currentColName]]@sample_designation)
    
    # Collapse table to figure out what is shared from infant's perspective -
    #   remove any ASVs that have a count of NA in the current infant.
    F1Table <- F1Table[!(is.na(F1Table[[F1Designation]])), ]
    
    # Add a column to indicate if each ASV (row) in the collapsed table is 
    #   shared with mom or not. If the ASV has a value of NA in mother column
    #   then it cannot be shared between current infant and mom.
    
    F1Table$shared <- ifelse(!(is.na(F1Table$"Manure-Manure")), yes = TRUE, no = FALSE)
    F1Table$sample <- F1Identifier

    # Prepare to add this table to a bigger list, but don't add until end!
    setnames(F1Table, old = F1Designation, new = "read_count")
      
    # For this infant, calculate percent of reads that are from ASVs that
    #   are shared with mom, and not shared with mom (keep this table).
    F1Sum <- sum(F1Table$read_count)
    F1SharedSum <- sum(F1Table$read_count[F1Table$shared == TRUE])
    F1NotSharedSum <- sum(F1Table$read_count[F1Table$shared == FALSE])
      
    # Add the infantTable (the one with ASVs as rows) to larger data frame
    F1TablePercentages <- F1Table %>%
      mutate("F1_percent" = read_count/F1Sum) %>%
      dplyr::select(-c("Manure-Manure"))
    
    # Add infant sample name to table
    #infantTablePercentages <- cbind(infantTablePercentages,
                                    #"infant_sample" = currentInfantName)
    
    allF1ASVTables <- rbind(allF1ASVTables, F1TablePercentages)
  }
}

# Infant Perspective Richness

# Generate a plot that is similar to the one above, but by richness instead
#   of relative abundance of reads.

# Group allInfantASVTables by sample and shared (T/F) then count

F1PerspectiveRichness <- allF1ASVTables
array <- unlist(strsplit(F1PerspectiveRichness$sample, "[/]"))
date <- array[c(TRUE, FALSE)]
F1PerspectiveRichness$date <- date

F1PerspectiveRichness.1 <- F1PerspectiveRichness %>%
  group_by(date, sample, shared) %>%
  summarise(number_of_ASVs = n())
  
lst <- strsplit(F1PerspectiveRichness.1$sample,'/')
v2 <- lapply(lst, `[`, 2)
v2 = substr(v2,1,4)
F1PerspectiveRichness.1$sample.type <- v2

F1SharedOTUs <- F1PerspectiveRichness.1
F1SharedOTUs <- filter(F1SharedOTUs, shared == TRUE)
F1SharedOTUs <- subset(F1SharedOTUs, select = c(date,sample,number_of_ASVs, sample.type))
F1SharedOTUs$sample.type <- factor(F1SharedOTUs$sample.type, levels=c("Endo", "Ecto"))
kruskal.test(number_of_ASVs~as.factor(sample.type),data=F1SharedOTUs) #Significant (p < 0.0001)
means <- aggregate(number_of_ASVs ~ sample.type, data = F1SharedOTUs, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(number_of_ASVs ~ sample.type, data = F1SharedOTUs, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="sample.type")
summary$sample.type <- factor(summary$sample.type, levels=c("Endo", "Ecto"))

#### Fig. S7 (right, endo)    
F1SharedOTUs$date <- factor(F1SharedOTUs$date, levels=c("7.8.2021","7.15.2021","7.22.2021","7.29.2021","8.5.2021","8.12.2021","8.19.2021","8.26.2021","9.9.2021","9.12.2021","9.15.2021"))
F1SharedOTUs.Endo <- F1SharedOTUs[F1SharedOTUs$sample.type=="Endo",]
F1SharedOTUs.Ecto <- F1SharedOTUs[F1SharedOTUs$sample.type=="Ecto",]
kruskal.test(number_of_ASVs~as.factor(date),data=F1SharedOTUs.Endo) #NS (p = 0.5275)
kruskal.test(number_of_ASVs~as.factor(date),data=F1SharedOTUs.Ecto) #NS (p = 0.0781)

means <- aggregate(number_of_ASVs ~ date, data = F1SharedOTUs.Endo, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(number_of_ASVs ~ date, data = F1SharedOTUs.Endo, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="date")
summary$date <- factor(summary$date, levels=c("7.8.2021","7.15.2021","7.22.2021","7.29.2021","8.5.2021","8.12.2021","8.19.2021","8.26.2021","9.9.2021","9.12.2021","9.15.2021"))
ggplot(summary, aes(x = date, y = number_of_ASVs.x)) +
    geom_errorbar(aes(ymin=number_of_ASVs.x-number_of_ASVs.y, ymax=number_of_ASVs.x+number_of_ASVs.y), width=.1) + #Fig. S7 (left)
    geom_point()

#### Fig. S7 (right, ecto)  
means <- aggregate(number_of_ASVs ~ date, data = F1SharedOTUs.Ecto, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(number_of_ASVs ~ date, data = F1SharedOTUs.Ecto, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="date")
summary$date <- factor(summary$date, levels=c("7.8.2021","7.15.2021","7.22.2021","7.29.2021","8.5.2021","8.12.2021","8.19.2021","8.26.2021","9.9.2021","9.12.2021","9.15.2021"))
ggplot(summary, aes(x = date, y = number_of_ASVs.x)) +
    geom_errorbar(aes(ymin=number_of_ASVs.x-number_of_ASVs.y, ymax=number_of_ASVs.x+number_of_ASVs.y), width=.1) + #Fig. S7 (left)
    geom_point()

##

F1PerspectiveRA <- allF1ASVTables
array <- unlist(strsplit(F1PerspectiveRA$sample, "[/]"))
date <- array[c(TRUE, FALSE)]
F1PerspectiveRA$date <- date

F1PerspectiveRA.1 <- F1PerspectiveRA %>%
  group_by(date, sample, shared) %>%
  summarise(tot_shared_reads = sum(read_count))
  
F1PerspectiveRA.2 <- F1PerspectiveRA.1 %>%
  group_by(sample) %>%
  summarise(tot_reads = sum(tot_shared_reads))

F1PerspectiveRA.3 <- merge(x = F1PerspectiveRA.1,y = F1PerspectiveRA.2,by.x = "sample",by.y = "sample", all = T)

F1PerspectiveRA.3$percent <- F1PerspectiveRA.3$tot_shared_reads/F1PerspectiveRA.3$tot_reads

write.csv(F1PerspectiveRA.3, 
            file = "ManurePerspectiveRAFlies.csv",
            quote = FALSE, row.names = FALSE)

#### Fig. S9 (top right, endo)
                 
data=read.csv("DCCManurePerspectiveRAFliesFormatted.csv")
data.endo <- data[data$fly.source=="Endo",]
data_summary <- aggregate(percent_shared ~ as.factor(date), data.endo,       # Create summary data
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("7.8.2021","7.15.2021","7.22.2021","7.29.2021","8.5.2021","8.12.2021","8.19.2021","8.26.2021","9.9.2021","9.12.2021","9.15.2021"))
ggplot(data_summary) + #Fig. S9
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)

#### Fig. S9 (bottom right, ecto)
                          
data.ecto <- data[data$fly.source=="Ecto",]
data_summary <- aggregate(percent_shared ~ as.factor(date), data.ecto,       # Create summary data
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("7.8.2021","7.15.2021","7.22.2021","7.29.2021","8.5.2021","8.12.2021","8.19.2021","8.26.2021","9.9.2021","9.12.2021","9.15.2021"))
ggplot(data_summary) + #Fig. S9
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)

### Fig. S11

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
sample_data(ps)$merge_factor <- paste(sample_data(ps)$sample.type,sample_data(ps)$sampling.date, sep = "-")
ps.arlington <- subset_samples(ps, location=="Arlington")
ps.Merged <- merge_samples(ps.arlington,"merge_factor")
row_names <- row.names(sample_data(ps.Merged))
sample_data(ps.Merged)$sample_id <- row_names
lst <- strsplit(row_names,'-')
v1 <- lapply(lst, `[`, 1)
v2 <- lapply(lst, `[`, 2)
sample_data(ps.Merged)$sample.type <- v1
sample_data(ps.Merged)$sampling.date <- v2

setClass("Family",
         representation = representation(family_group = "character",
                                         member_list = "list",
                                         family_sum = "numeric"),
         prototype = prototype(family_group = NA_character_,
                               member_list = list(),
                               family_sum = NA_real_))

setClass("FamilyMember", 
         representation = representation(family_group = "character",
                                         name = "character",
                                         sample_designation = "character",
                                         sample_sum = "numeric",
                                         ASV_counts = "data.frame"),
         prototype = prototype(family_group = NA_character_,
                               name = NA_character_,
                               sample_designation = NA_character_,
                               sample_sum = NA_real_))

CountASVs <- function(sampleName, abundanceTable) {
  sampleAbundance <- abundanceTable %>%
    dplyr::select(sampleName) %>%
    rownames_to_column(var = "ASV")
  sampleAbundanceDF <- filter(sampleAbundance, sampleAbundance[, sampleName] > 0)
  return(sampleAbundanceDF)
}

GetSampleSum <- function(sampleName, physeqObject) {
  named_sum <- sample_sums(physeqObject)[sampleName]
  numeric_sum <- unname(named_sum)
}

# Pull out sample data
sampleData <- data.frame(sample_data(ps.Merged))
sampleData$sample.type <- unlist(sampleData$sample.type)
sampleData$sampling.date <- unlist(sampleData$sampling.date)

# Pull out ASV count table
abdTable <- as.matrix(otu_table(ps.Merged))
if (!taxa_are_rows(abdTable)) {
  abdTable <- t(abdTable)
}
abdTable <- as.data.frame(abdTable)

# Populate Classes

# This list will hold the Family objects once completed
familyClassList <- list()

# Loop over the rows of the sample data data frame, building the objects by
#   familiy group number.  But we only want to build it once, so when we run
#   across the family group number again we don't want to re-build it.

for (i in 1:nrow(sampleData)) {
  # Get the current family group from the sample data table
  currentFamilyGroup <- sampleData[i, "sampling.date"]
  
  # Get all the current family names in familyClassList to check if 
  #   currentFamilyGroup has already been built:
  familyNames <- names(familyClassList)
  
  # If the current family group has not been built, build it:
  if (!(currentFamilyGroup %in% familyNames)){
    # Populate the family_group slot:
    # Create a list of sample names for all the samples in the current family group
    #   by subsetting the sample data df and extracting the sample names as a char
    #   vector.  Then, add to the memberList list that will be looped over to
    #   build each FamilyMember object for the current Family.
    members <- sampleData %>%
      filter(sampling.date == currentFamilyGroup) %>%
      dplyr::select(all_of("sample_id"))
    memberList <- as.list(members$sample_id)
    
    # Create an empty list to hold FamilyMember objects which will be added to the 
    #   Family object in the end.
    currentFamilyMembersList <- list()
    
    # Loop over the sample names in memberList and build a new FamilyMember object for
    #   each sample:
    for (sample in memberList) {
      newMember <- new("FamilyMember",
                       family_group = currentFamilyGroup,
                       name = sample, 
                       sample_designation = sampleData$sample.type[sampleData$sample_id == sample],
                       sample_sum = GetSampleSum(sample, ps.Merged),
                       ASV_counts = CountASVs(sample, abdTable))
    # Add this new FamilyMember object to currentFamilyMemberList:
    currentFamilyMembersList[[sample]] <- newMember
    }
    
    # Now that the member list is complete, get a total read sum to be put into the 
    #  Family object (sum of each family member's total read count).
    totalSum <- 0
    for (member in currentFamilyMembersList) {
      sampleSum <- member@sample_sum
      totalSum <- sum(totalSum, sampleSum)
    }
    
    # We now have all the components needed to create new Family object.
    newFamily <- new("Family", 
                     family_group = currentFamilyGroup,
                     member_list = currentFamilyMembersList,
                     family_sum = totalSum)
    # Add to list
    familyClassList[[currentFamilyGroup]] <- newFamily
  }
}

# Select "Valid" Families

#Families including a mother and at least one infant.

# A valid Family must have at least one mother and one infant:
# So a Family object should have a family_member list equal to length 3 or
# equal to length 2 where one of the sample_designation values is "mother"

# Create a list to hold the valid Family objects:
validFamilies <- familyClassList

# Infant Perspective 

# This is for the "infant perspective" plot.

# Goal: generate a data frame holding all infants and all their ASVs (as rows). For each
#   ASV in the given infant, count the number of reads and
#   the percentage of the total reads in the infant that the ASV accounts for
#   ("infant_percent" column.)

allF1ASVTables <- data.frame()

for (k in 1:length(validFamilies)) {
  family <- validFamilies[[k]]
  familyMembers <- family@member_list

  # Put each FamilyMember object in its own variable
  currentManure <- NULL
  currentEndo <- NULL
  currentEcto <- NULL

  for (i in 1:length(familyMembers)) {
    if (familyMembers[[i]]@sample_designation == "Manure") {
      currentManure <- familyMembers[[i]]
    } else if (familyMembers[[i]]@sample_designation == "Endo") {
      currentEndo <- familyMembers[[i]]
    } else if (familyMembers[[i]]@sample_designation == "Ecto") {
      currentEcto <- familyMembers[[i]]
    } 
  }

  # Put infants in a list to loop over:
  currentF1List <- list(currentManure, currentEndo, currentEcto)
  
  # Remove any NULL elements in the list (in case there is only 1 infant)
  currentF1List <- plyr::compact(currentF1List)

  # Loop over the infant list and populate the "infant perspective" data frame:
  #   allInfantASVTables - for each infant samples, holds the ASV, 
  #     its percentage of the total reads in the infant, and whether it's shared with mother
  #   (Later) to get the total percent shared and not shared for each infant (for plotting)
  #     allInfantASVTables %>% group_by(sample, shared) %>% summarise(percent = sum(infants))
  
  for (i in 1:length(currentF1List)) {
    
    current_F1 <- currentF1List[[i]]
    #currentInfantName <- current_infant@name
    
    # build infant's identifier (FamGroup.1 or FamGroup.2)
    F1Identifier <- paste(current_F1@sample_designation, current_F1@family_group, sep = "-")
    F1Table <- current_F1@ASV_counts
    
    # Replace column names with the samples' sample_designations for easier
    #   plotting downstream:
#    for (j in 1:ncol(F1Table)) {
#      currentColName <- names(F1Table)[j]
#      if (currentColName != "ASV") {
#        setnames(F1Table,
#                 old = currentColName,
#                 new = familyMembers[[currentColName]]@sample_designation)
#      }
#    }
                
    # Prepare to add this table to a bigger list, but don't add until end!
    setnames(F1Table, old = F1Identifier, new = "read_count", skip_absent = TRUE)
      
    # For this infant, calculate percent of reads that are from ASVs that
    #   are shared with mom, and not shared with mom (keep this table).
    F1Sum <- sum(F1Table$read_count)
      
    # Add the infantTable (the one with ASVs as rows) to larger data frame
    F1TablePercentages <- F1Table %>%
      mutate("F1_percent" = read_count/F1Sum)
    
    # Add infant sample name to table
    F1TablePercentages <- cbind(F1TablePercentages,
                                    "sample" = F1Identifier)
    
    allF1ASVTables <- rbind(allF1ASVTables, F1TablePercentages)
  }
}

allF1ASVTables <- allF1ASVTables[allF1ASVTables$ASV %in% c("15d888a8d8e6c158974710b615dcbcdf",
"ff27974c65664de75e6421a62a0448a3",
"b40dfe8590cdabab904a9931210d2805",
"b77e07971c99cc0a36c99bd7cad83c16",
"5dbd59fa8200705770f4d068d6572dc3",
"bcaaa3d15ff438cfa14dbec5d86b0eb8",
"e0597b0224419e26dd82616c75f8c1d2",
"945bb343f8dc3dd656b9ca15c3c0f09d",
"dcfd8312bf9f39a54aab3707c3b140df",
"350bce411d627b38ef18aeaebd7d4f69",
"91f2efe0637f2f2062c1252c32a8e4eb",
"6207db2806b88c7f75541b8de1c42f50",
"66402d8adb61bb82bd8d7fb58d832c10",
"8da7704aad448995f55464ee650a4452",
"3289d95c4a81c516faa9d263f55202d2",
"aaa9219facfcb711b04e0d3c190b571a",
"d56a560e43e7ab41839a1c8a10f65d36"), ]

allF1ASVTables$sample <- factor(allF1ASVTables$sample,levels = c("Endo-7.9.2021", "Manure-7.9.2021",
	"Endo-7.16.2021", "Ecto-7.16.2021", "Manure-7.16.2021",
	"Endo-7.23.2021", "Ecto-7.23.2021", "Manure-7.23.2021",
	"Endo-7.30.2021", "Ecto-7.30.2021", "Manure-7.30.2021",
	"Endo-8.6.2021", "Ecto-8.6.2021", "Manure-8.6.2021",
	"Endo-8.13.2021", "Ecto-8.13.2021", "Manure-8.13.2021",
	"Endo-8.20.2021", "Ecto-8.20.2021", "Manure-8.20.2021",
	"Endo-8.27.2021", "Ecto-8.27.2021", "Manure-8.27.2021",
	"Endo-9.10.2021", "Ecto-9.10.2021", "Manure-9.10.2021"))

ASVTax <- as.data.frame(tax_table(ps.Merged))
ASVIDs <- row.names(ASVTax)
ASVTax$ASV <- ASVIDs
ASVTax <- ASVTax[ASVTax$ASV %in% c("15d888a8d8e6c158974710b615dcbcdf",
"ff27974c65664de75e6421a62a0448a3",
"b40dfe8590cdabab904a9931210d2805",
"b77e07971c99cc0a36c99bd7cad83c16",
"5dbd59fa8200705770f4d068d6572dc3",
"bcaaa3d15ff438cfa14dbec5d86b0eb8",
"e0597b0224419e26dd82616c75f8c1d2",
"945bb343f8dc3dd656b9ca15c3c0f09d",
"dcfd8312bf9f39a54aab3707c3b140df",
"350bce411d627b38ef18aeaebd7d4f69",
"91f2efe0637f2f2062c1252c32a8e4eb",
"6207db2806b88c7f75541b8de1c42f50",
"66402d8adb61bb82bd8d7fb58d832c10",
"8da7704aad448995f55464ee650a4452",
"3289d95c4a81c516faa9d263f55202d2",
"aaa9219facfcb711b04e0d3c190b571a",
"d56a560e43e7ab41839a1c8a10f65d36"), ]

#	Domain	Phylum	Class	Order	Family	Genus
#3289d95c4a81c516faa9d263f55202d2	Bacteria	Proteobacteria	Gammaproteobacteria	Pseudomonadales	Moraxellaceae	Acinetobacter
#350bce411d627b38ef18aeaebd7d4f69	Bacteria	Proteobacteria	Gammaproteobacteria	Pseudomonadales	Moraxellaceae	Acinetobacter
#6207db2806b88c7f75541b8de1c42f50	Bacteria	Proteobacteria	Gammaproteobacteria	Enterobacterales	Morganellaceae	Providencia
#91f2efe0637f2f2062c1252c32a8e4eb	Bacteria	Proteobacteria	Gammaproteobacteria	Enterobacterales	Enterobacteriaceae	Escherichia-Shigella
#5dbd59fa8200705770f4d068d6572dc3	Bacteria	Proteobacteria	Gammaproteobacteria	Enterobacterales	Erwiniaceae	Pantoea

#66402d8adb61bb82bd8d7fb58d832c10	Bacteria	Firmicutes	Bacilli	Lactobacillales	Enterococcaceae	Enterococcus
#d56a560e43e7ab41839a1c8a10f65d36	Bacteria	Firmicutes	Bacilli	Lactobacillales	Enterococcaceae	Enterococcus
#8da7704aad448995f55464ee650a4452	Bacteria	Firmicutes	Bacilli	Staphylococcales	Staphylococcaceae	Jeotgalicoccus
#ff27974c65664de75e6421a62a0448a3	Bacteria	Firmicutes	Bacilli	Staphylococcales	Staphylococcaceae	Staphylococcus
#15d888a8d8e6c158974710b615dcbcdf	Bacteria	Firmicutes	Bacilli	Staphylococcales	Staphylococcaceae	Staphylococcus
#945bb343f8dc3dd656b9ca15c3c0f09d	Bacteria	Firmicutes	Bacilli	Staphylococcales	Staphylococcaceae	Staphylococcus
#e0597b0224419e26dd82616c75f8c1d2	Bacteria	Firmicutes	Bacilli	Lactobacillales	Vagococcaceae	Vagococcus

#b40dfe8590cdabab904a9931210d2805	Bacteria	Actinobacteriota	Actinobacteria	Micrococcales	Brevibacteriaceae	Brevibacterium
#bcaaa3d15ff438cfa14dbec5d86b0eb8	Bacteria	Actinobacteriota	Actinobacteria	Micrococcales	Dermabacteraceae	Brachybacterium
#aaa9219facfcb711b04e0d3c190b571a	Bacteria	Actinobacteriota	Actinobacteria	Corynebacteriales	Dietziaceae	Dietzia
#b77e07971c99cc0a36c99bd7cad83c16	Bacteria	Actinobacteriota	Actinobacteria	Micrococcales	Intrasporangiaceae	Unidentified
#dcfd8312bf9f39a54aab3707c3b140df	Bacteria	Actinobacteriota	Actinobacteria	Micrococcales	Micrococcaceae	Glutamicibacter

#6721
#6741
#6852
#6844
#6837

#7636
#7651
#7469
#7480
#7491
#7497
#7614

#4904
#5078
#4754
#5106
#4966


allF1ASVTables$ASV <- factor(allF1ASVTables$ASV,levels = c("3289d95c4a81c516faa9d263f55202d2",
"350bce411d627b38ef18aeaebd7d4f69",
"6207db2806b88c7f75541b8de1c42f50",
"91f2efe0637f2f2062c1252c32a8e4eb",
"5dbd59fa8200705770f4d068d6572dc3",
"66402d8adb61bb82bd8d7fb58d832c10",
"d56a560e43e7ab41839a1c8a10f65d36",
"8da7704aad448995f55464ee650a4452",
"ff27974c65664de75e6421a62a0448a3",
"15d888a8d8e6c158974710b615dcbcdf",
"945bb343f8dc3dd656b9ca15c3c0f09d",
"e0597b0224419e26dd82616c75f8c1d2",
"b40dfe8590cdabab904a9931210d2805",
"bcaaa3d15ff438cfa14dbec5d86b0eb8",
"aaa9219facfcb711b04e0d3c190b571a",
"b77e07971c99cc0a36c99bd7cad83c16",
"dcfd8312bf9f39a54aab3707c3b140df"))

ASVNames <- c("3289d95c4a81c516faa9d263f55202d2",
"350bce411d627b38ef18aeaebd7d4f69",
"6207db2806b88c7f75541b8de1c42f50",
"91f2efe0637f2f2062c1252c32a8e4eb",
"5dbd59fa8200705770f4d068d6572dc3",
"66402d8adb61bb82bd8d7fb58d832c10",
"d56a560e43e7ab41839a1c8a10f65d36",
"8da7704aad448995f55464ee650a4452",
"ff27974c65664de75e6421a62a0448a3",
"15d888a8d8e6c158974710b615dcbcdf",
"945bb343f8dc3dd656b9ca15c3c0f09d",
"e0597b0224419e26dd82616c75f8c1d2",
"b40dfe8590cdabab904a9931210d2805",
"bcaaa3d15ff438cfa14dbec5d86b0eb8",
"aaa9219facfcb711b04e0d3c190b571a",
"b77e07971c99cc0a36c99bd7cad83c16",
"dcfd8312bf9f39a54aab3707c3b140df")

#nb.cols <- length(ASVNames)
#myColors <- colorRampPalette(brewer.pal(n = 8, name = "Set2"))(nb.cols)
myColors <- c("#edf3f7",
"#d2e2ef",
"#9fc1dc",
"#4292c6",
"#2171b5",
"#fee6ce",
"#fdae6b",
"#fd8d3c",
"#f16913",
"#d94801",
"#a63603",
"#7f2704",
"#f1efe1",
"#e9e9e7",
"#d4d4d4",
"#bababa",
"#878787")
names(myColors) <- ASVNames

F1PerspectivePlot <- ggplot(allF1ASVTables,
                                aes(x = sample, y = F1_percent, 
                                    fill = ASV)) +
  geom_bar(stat = "identity") +
scale_fill_manual(name = "ASV", values = myColors)

F1PerspectivePlot #Fig. S11
