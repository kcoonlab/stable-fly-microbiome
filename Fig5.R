set.seed(123456)

# Load required packages

# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")

library(phyloseq)
library(dplyr)
library(tidyverse)
library(data.table)

#### Fig 5A (left, Arlington samples)

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
    dplyr::select(all_of(sampleName)) %>%
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

# Loop over the rows of the sample data data frame, building the objects by "familiy" group number.  But we only want to build it once, so when we run across the family group number again we don't want to re-build it.

for (i in 1:140) {
  # Get the current family group from the sample data table
  currentFamilyGroup <- sampleData[i, "sampling.date"]
  
  # Get all the current family names in familyClassList to check if currentFamilyGroup has already been built:
  familyNames <- names(familyClassList)
  
  # If the current family group has not been built, build it:
  if (!(currentFamilyGroup %in% familyNames)){
    # Populate the family_group slot:
    # Create a list of sample names for all the samples in the current family group by subsetting the sample data df and extracting the sample names as a char vector.  Then, add to the memberList list that will be looped over to build each FamilyMember object for the current Family.
    members <- sampleData %>%
      filter(sampling.date == currentFamilyGroup) %>%
      dplyr::select(all_of("sample_id"))
    memberList <- as.list(members$sample_id)
    
    # Create an empty list to hold FamilyMember objects which will be added to the Family object in the end.
    currentFamilyMembersList <- list()
    
    # Loop over the sample names in memberList and build a new FamilyMember object for each sample:
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
    
    # Now that the member list is complete, get a total read sum to be put into the Family object (sum of each family member's total read count).
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

# A valid Family must have at least one manure sample and one fly sample:
# So a Family object should have a family_member list equal to length 3 or equal to length 2 where one of the sample_designation values is "Manure"

# Create a list to hold the valid Family objects:
validFamilies <- familyClassList

# Goal: generate a data frame holding all fly samples and all their ASVs (as rows). 
# For each ASV in the given fly sample, count the number of reads and the percentage of the total reads in the fly sample that the ASV accounts for

allFlyASVTables <- data.frame()

for (k in 1:length(validFamilies)) {
  family <- validFamilies[[k]]
  familyMembers <- family@member_list

  # Put each FamilyMember object in its own variable
  currentFly <- NULL

  for (i in 1:length(familyMembers)) {
      currentFly <- c(currentFly,familyMembers[[i]])
    	}

  # Put fly samples in a list to loop over:
  currentFlyList <- list(currentFly)
  
  # Remove any NULL elements in the list (in case there is only 1 fly sample)
  currentFlyList <- plyr::compact(currentFlyList)

  # Loop over the fly sample list and populate the "fly perspective" data frame:
  #   allFlyASVTables - for each fly sample, holds the ASV, its percentage of the total reads in the sample, and whether it's shared with manure
  
  for (i in 1:length(currentFlyList[[1]])) {
  
    current_Fly <- currentFlyList[[1]][[i]]
    
    # build fly sample's identifier (FamGroup.1 or FamGroup.2)
    FlyDesignation <- current_Fly@sample_designation
    FlyIdentifier <- paste(family@family_group, FlyDesignation, sep = "/")
	
	Manure_ASV_counts = CountASVs("Manure-Manure", abdTable)

    # merge current fly sample's ASV counts with manure's
    FlyTable <- merge(current_Fly@ASV_counts,
                         Manure_ASV_counts,
                         by = "ASV", all = TRUE)
    
    # Replace column names with the samples' sample_designations for easier plotting downstream:
    
      currentColName <- names(FlyTable)[2]
      setnames(FlyTable, old = currentColName, new = familyMembers[[currentColName]]@sample_designation)
    
    # Collapse table to figure out what is shared from fly sample's perspective - remove any ASVs that have a count of NA in the current sample
    FlyTable <- FlyTable[!(is.na(FlyTable[[FlyDesignation]])), ]
    
    # Add a column to indicate if each ASV (row) in the collapsed table is shared with manure or not. If the ASV has a value of NA in manure column then it cannot be shared between current fly sample and manure
    
    FlyTable$shared <- ifelse(!(is.na(FlyTable$"Manure-Manure")), yes = TRUE, no = FALSE)
    FlyTable$sample <- FlyIdentifier

    # Prepare to add this table to a bigger list, but don't add until end!
    setnames(FlyTable, old = FlyDesignation, new = "read_count")
      
    # For this fly sample, calculate percent of reads that are from ASVs that are shared with manure, and not shared with manure (keep this table).
    FlySum <- sum(FlyTable$read_count)
    FlySharedSum <- sum(FlyTable$read_count[FlyTable$shared == TRUE])
    FlyNotSharedSum <- sum(FlyTable$read_count[FlyTable$shared == FALSE])
      
    # Add the flyTable (the one with ASVs as rows) to larger data frame
    FlyTablePercentages <- FlyTable %>%
      mutate("Fly_percent" = read_count/FlySum) %>%
      dplyr::select(-c("Manure-Manure"))
        
    allFlyASVTables <- rbind(allFlyASVTables, FlyTablePercentages)
  }
}

# Generate a plot that is similar to the one above, but by richness instead of relative abundance of reads.

# Group allFlyASVTables by sample and shared (T/F) then count

FlyPerspectiveRichness <- allFlyASVTables
array <- unlist(strsplit(FlyPerspectiveRichness$sample, "[/]"))
date <- array[c(TRUE, FALSE)]
FlyPerspectiveRichness$date <- date

FlyPerspectiveRichness.1 <- FlyPerspectiveRichness %>%
  group_by(date, sample, shared) %>%
  summarise(number_of_ASVs = n())
  
lst <- strsplit(FlyPerspectiveRichness.1$sample,'/')
v2 <- lapply(lst, `[`, 2)
v2 = substr(v2,1,4)
FlyPerspectiveRichness.1$sample.type <- v2

FlySharedOTUs <- FlyPerspectiveRichness.1
FlySharedOTUs <- filter(FlySharedOTUs, shared == TRUE)
FlySharedOTUs <- subset(FlySharedOTUs, select = c(date,sample,number_of_ASVs, sample.type))
FlySharedOTUs$sample.type <- factor(FlySharedOTUs$sample.type, levels=c("Endo", "Ecto"))
FlySharedOTUs$trap <- c("B2-West","B1-East","B2-East","B2-West","B2-West","B2-West","B1-East","B1-East","B1-East","B1-East","B1-West","B1-West","B1-West","B1-West","B2-East","B2-East","B2-East","B1-East","B1-West","B2-East","B1-East","B1-East","B1-East","B1-East","B1-West","B1-West","B1-West","B1-West","B1-West","B2-East","B2-East","B2-East","B2-West","B2-West","B2-West","B2-West","B1-East","B1-West","B2-East","B2-West","B1-East","B1-West","B1-West","B2-East","B2-East","B2-West","B2-West","B1-East","B1-East","B1-West","B1-West","B1-West","B1-West","B2-East","B2-East","B2-East","B2-East","B2-West","B2-West","B1-East","B1-West","B2-East","B2-West","B1-East","B1-East","B1-East","B1-East","B1-East","B1-West","B1-West","B2-East","B2-East","B2-West","B2-West","B2-West","B2-West","B1-East","B1-West","B2-East","B2-West","B1-East","B1-East","B1-West","B1-West","B2-East","B2-East","B2-East","B2-West","B1-East","B1-West","B2-East","B2-West","B1-East","B1-East","B1-East","B1-East","B1-East","B1-West","B2-East","B2-East","B2-East","B2-East","B2-West","B2-West","B2-West","B1-East","B1-West","B2-East","B2-West","B1-East","B1-East","B1-East","B1-East","B1-East","B1-West","B1-West","B1-West","B1-West","B1-West","B2-East","B2-East","B2-East","B2-East","B2-East","B2-West","B2-West","B2-West","B2-West","B1-East","B2-East","B1-East","B1-East","B1-East","B1-East","B1-East","B2-East","B2-East","B2-East","B2-East","B2-East")
FlySharedOTUs$trap <- factor(FlySharedOTUs$trap, levels=c("B1-East", "B1-West", "B2-East", "B2-West"))
kruskal.test(number_of_ASVs~as.factor(sample.type),data=FlySharedOTUs) #Significant (p < 0.0001)
means <- aggregate(number_of_ASVs ~ sample.type, data = FlySharedOTUs, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(number_of_ASVs ~ sample.type, data = FlySharedOTUs, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="sample.type")
summary$sample.type <- factor(summary$sample.type, levels=c("Endo", "Ecto"))

ggplot(summary, aes(x = sample.type, y = number_of_ASVs.x)) + #Fig. 5A
    geom_errorbar(aes(ymin=number_of_ASVs.x-number_of_ASVs.y, ymax=number_of_ASVs.x+number_of_ASVs.y), width=.1) +
    geom_point()
  
#### Fig 5B (left, Arlington samples)

FlyPerspectiveRA <- allFlyASVTables
array <- unlist(strsplit(FlyPerspectiveRA$sample, "[/]"))
date <- array[c(TRUE, FALSE)]
FlyPerspectiveRA$date <- date

FlyPerspectiveRA.1 <- FlyPerspectiveRA %>%
  group_by(date, sample, shared) %>%
  summarise(tot_shared_reads = sum(read_count))
  
FlyPerspectiveRA.2 <- FlyPerspectiveRA.1 %>%
  group_by(sample) %>%
  summarise(tot_reads = sum(tot_shared_reads))

FlyPerspectiveRA.3 <- merge(x = FlyPerspectiveRA.1,y = FlyPerspectiveRA.2,by.x = "sample",by.y = "sample", all = T)

FlyPerspectiveRA.3$percent <- FlyPerspectiveRA.3$tot_shared_reads/FlyPerspectiveRA.3$tot_reads

#write.csv(F1PerspectiveRA.3, file = "ManurePerspectiveRAFlies.csv", quote = FALSE, row.names = FALSE)

data=read.csv("ManurePerspectiveRAFliesFormatted.csv")
kruskal.test(percent_shared~as.factor(fly.source),data=data) #Significant (p = 0.04489)
kruskal.test(percent_shared~as.factor(date),data=data) #Significant (p = 0.01791)
pairwise.wilcox.test(data$percent_shared, as.factor(data$date),
                 p.adjust.method = "BH", exact = FALSE) #Only significant: 7.9.2021 vs. 8.13.2021

data_summary <- aggregate(percent_shared ~ as.factor(fly.source), data,       # Create summary data
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("Endo", "Ecto"))
ggplot(data_summary) + #Fig. 5B
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, linewidth=1.3) +
    ylim(0,1)

#### Fig 5C (left, Arlington samples)

ps <- readRDS("ps_FieldWork2021_AJS_Final.rds")
sample_data(ps)$merge_factor <- c("Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","ArlM1-7.9.21","ArlM10-7.16.21","ArlM11-7.23.21","ArlM12-7.23.21","ArlM13-7.23.21","ArlM14-7.23.21","ArlM15-7.23.21","ArlM16-7.30.21","ArlM17-7.30.21","ArlM18-7.30.21","ArlM19-7.30.21","ArlM2-7.9.21","ArlM20-7.30.21","ArlM21-7.30.21","ArlM22-8.6.21","ArlM23-8.6.21","ArlM24-8.6.21","ArlM25-8.6.21","ArlM26-8.6.21","ArlM27-8.6.21","ArlM28-8.13.21","ArlM29-8.13.21","ArlM3-7.9.21","ArlM30-8.13.21","ArlM31-8.13.21","ArlM32-8.13.21","ArlM33-8.13.21","ArlM34-8.20.21","ArlM35-8.20.21","ArlM36-8.20.21","ArlM37-8.20.21","ArlM38-8.20.21","ArlM39-8.20.21","ArlM4-7.9.21","ArlM40-8.27.21","ArlM41-8.27.21","ArlM42-8.27.21","ArlM43-8.27.21","ArlM44-8.27.21","ArlM45-8.27.21","ArlM46-9.10.21","ArlM47-9.10.21","ArlM48-9.10.21","ArlM49-9.10.21","ArlM5-7.9.21","ArlM50-9.10.21","ArlM51-9.10.21","ArlM6-7.16.21","ArlM7-7.16.21","ArlM8-7.16.21","ArlM9-7.16.21","DCCM100-9.15.21","DCCM101-9.15.21","DCCM102-7.8.21","DCCM103-7.8.21","DCCM104-7.8.21","DCCM105-7.8.21","DCCM106-7.8.21","DCCM52-7.15.21","DCCM53-7.15.21","DCCM54-7.15.21","DCCM55-7.15.21","DCCM56-7.15.21","DCCM57-7.22.21","DCCM58-7.22.21","DCCM59-7.22.21","DCCM60-7.22.21","DCCM61-7.22.21","DCCM62-7.29.21","DCCM63-7.29.21","DCCM64-7.29.21","DCCM65-7.29.21","DCCM66-7.29.21","DCCM67-8.5.21","DCCM68-8.5.21","DCCM69-8.5.21","DCCM70-8.5.21","DCCM71-8.5.21","DCCM72-8.12.21","DCCM73-8.12.21","DCCM74-8.12.21","DCCM75-8.12.21","DCCM76-8.12.21","DCCM77-8.19.21","DCCM78-8.19.21","DCCM79-8.19.21","DCCM80-8.19.21","DCCM81-8.19.21","DCCM82-8.26.21","DCCM83-8.26.21","DCCM84-8.26.21","DCCM85-8.26.21","DCCM86-8.26.21","DCCM87-9.2.21","DCCM88-9.2.21","DCCM89-9.2.21","DCCM90-9.2.21","DCCM91-9.2.21","DCCM92-9.9.21","DCCM93-9.9.21","DCCM94-9.9.21","DCCM95-9.9.21","DCCM96-9.9.21","DCCM97-9.15.21","DCCM98-9.15.21","DCCM99-9.15.21","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly")
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

# Loop over the rows of the sample data data frame, building the objects by familiy group number.  But we only want to build it once, so when we run across the family group number again we don't want to re-build it.

for (i in 1:51) {
  # Get the current family group from the sample data table
  currentFamilyGroup <- sampleData[i, "sampling.date"]
  
  # Get all the current family names in familyClassList to check if 
  #   currentFamilyGroup has already been built:
  familyNames <- names(familyClassList)
  
  # If the current family group has not been built, build it:
  if (!(currentFamilyGroup %in% familyNames)){
    # Populate the family_group slot:
    # Create a list of sample names for all the samples in the current family group by subsetting the sample data df and extracting the sample names as a char vector.  Then, add to the memberList list that will be looped over to build each FamilyMember object for the current Family.
    members <- sampleData %>%
      filter(sampling.date == currentFamilyGroup) %>%
      dplyr::select(all_of("sample_id"))
    memberList <- as.list(members$sample_id)
    
    # Create an empty list to hold FamilyMember objects which will be added to the Family object in the end.
    currentFamilyMembersList <- list()
    
    # Loop over the sample names in memberList and build a new FamilyMember object for each sample:
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
    
    # Now that the member list is complete, get a total read sum to be put into the Family object (sum of each family member's total read count).
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

# A valid Family must have at least one manure sample and one fly sample:
# So a Family object should have a family_member list equal to length 3 or equal to length 2 where one of the sample_designation values is "Manure"

# Create a list to hold the valid Family objects:
validFamilies <- familyClassList

# Goal: generate a data frame holding all manure samples and all their ASVs (as rows). For each ASV in the given manure sample, count the number of reads and the percentage of the total reads in the manure sample that the ASV accounts for.

allManureASVTables <- data.frame()

for (k in 1:length(validFamilies)) {
  family <- validFamilies[[k]]
  familyMembers <- family@member_list

  # Put each FamilyMember object in its own variable
  currentManure <- NULL

  for (i in 1:length(familyMembers)) {
      currentManure <- c(currentManure,familyMembers[[i]])
    	}

  # Put manure samples in a list to loop over:
  currentManureList <- list(currentManure)
  
  # Remove any NULL elements in the list (in case there is only 1 fly sample)
  currentManureList <- plyr::compact(currentManureList)

  # Loop over the manure sample list and populate the "manure perspective" data frame:
  #   allManureASVTables - for each manure sample, holds the ASV, its percentage of the total reads in the manure sample, and whether it's shared with a fly sample
  
  for (i in 1:length(currentManureList[[1]])) {
  
    current_Manure <- currentManureList[[1]][[i]]
    
    # build manure sample's identifier (FamGroup.1 or FamGroup.2)
    ManureDesignation <- current_Manure@sample_designation
    ManureIdentifier <- paste(family@family_group, ManureDesignation, sep = "/")
	
	Fly_ASV_counts = CountASVs("Fly-Fly", abdTable)

    # merge current manure sample's ASV counts with flies'
    ManureTable <- merge(current_Manure@ASV_counts,
                         Fly_ASV_counts,
                         by = "ASV", all = TRUE)
    
    # Replace column names with the samples' sample_designations for easier plotting downstream:
    
      currentColName <- names(ManureTable)[2]
      setnames(ManureTable, old = currentColName, new = familyMembers[[currentColName]]@sample_designation)
    
    # Collapse table to figure out what is shared from manure sample's perspective - remove any ASVs that have a count of NA in the current manure sample
    ManureTable <- ManureTable[!(is.na(ManureTable[[ManureDesignation]])), ]
    
    # Add a column to indicate if each ASV (row) in the collapsed table is shared with a fly sample or not. If the ASV has a value of NA in fly column then it cannot be shared between current manure sample and flies.
    
    ManureTable$shared <- ifelse(!(is.na(ManureTable$"Fly-Fly")), yes = TRUE, no = FALSE)
    ManureTable$sample <- ManureIdentifier

    # Prepare to add this table to a bigger list, but don't add until end!
    setnames(ManureTable, old = ManureDesignation, new = "read_count")
      
    # For this manure sample, calculate percent of reads that are from ASVs that are shared with flies, and not shared with flies (keep this table).
    ManureSum <- sum(ManureTable$read_count)
    ManureSharedSum <- sum(ManureTable$read_count[ManureTable$shared == TRUE])
    ManureNotSharedSum <- sum(ManureTable$read_count[ManureTable$shared == FALSE])
      
    # Add the ManureTable (the one with ASVs as rows) to larger data frame
    ManureTablePercentages <- ManureTable %>%
      mutate("Manure_percent" = read_count/ManureSum) %>%
      dplyr::select(-c("Fly-Fly"))
    
    allManureASVTables <- rbind(allManureASVTables, ManureTablePercentages)
  }
}

ManurePerspectiveRA <- allManureASVTables
array <- unlist(strsplit(ManurePerspectiveRA$sample, "[/]"))
date <- array[c(TRUE, FALSE)]
ManurePerspectiveRA$date <- date

ManurePerspectiveRA.1 <- ManurePerspectiveRA %>%
  group_by(date, sample, shared) %>%
  summarise(tot_shared_reads = sum(read_count))
  
ManurePerspectiveRA.2 <- ManurePerspectiveRA.1 %>%
  group_by(sample) %>%
  summarise(tot_reads = sum(tot_shared_reads))

ManurePerspectiveRA.3 <- merge(x = ManurePerspectiveRA.1,y = ManurePerspectiveRA.2,by.x = "sample",by.y = "sample", all = T)

ManurePerspectiveRA.3$percent <- ManurePerspectiveRA.3$tot_shared_reads/ManurePerspectiveRA.3$tot_reads

#write.csv(ManurePerspectiveRA.3, file = "FlyPerspectiveRAManure.csv", quote = FALSE, row.names = FALSE)

data=read.csv("FlyPerspectiveRAManureFormatted.csv")
kruskal.test(percent_shared~as.factor(date),data=data) #NS (p = 0.9071)

mean <- mean(data$percent_shared)
ci <- 1.96*(sd(data$percent_shared) / sqrt(length(data$percent_shared)))                              
ggplot(data=NULL) + #Fig. 5C
    geom_bar( aes(x="", y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x="", ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)


#### Fig 5A (right, DCC samples)

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
    dplyr::select(all_of(sampleName)) %>%
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

# Loop over the rows of the sample data data frame, building the objects by "familiy" group number.  But we only want to build it once, so when we run across the family group number again we don't want to re-build it.

for (i in 1:18) {
  # Get the current family group from the sample data table
  currentFamilyGroup <- sampleData[i, "sampling.date"]
  
  # Get all the current family names in familyClassList to check if currentFamilyGroup has already been built:
  familyNames <- names(familyClassList)
  
  # If the current family group has not been built, build it:
  if (!(currentFamilyGroup %in% familyNames)){
    # Populate the family_group slot:
    # Create a list of sample names for all the samples in the current family group by subsetting the sample data df and extracting the sample names as a char vector.  Then, add to the memberList list that will be looped over to build each FamilyMember object for the current Family.
    members <- sampleData %>%
      filter(sampling.date == currentFamilyGroup) %>%
      dplyr::select(all_of("sample_id"))
    memberList <- as.list(members$sample_id)
    
    # Create an empty list to hold FamilyMember objects which will be added to the Family object in the end.
    currentFamilyMembersList <- list()
    
    # Loop over the sample names in memberList and build a new FamilyMember object for each sample:
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
    
    # Now that the member list is complete, get a total read sum to be put into the Family object (sum of each family member's total read count).
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

# A valid Family must have at least one manure sample and one fly sample:
# So a Family object should have a family_member list equal to length 3 or equal to length 2 where one of the sample_designation values is "Manure"

# Create a list to hold the valid Family objects:
validFamilies <- familyClassList

# Goal: generate a data frame holding all fly samples and all their ASVs (as rows). 
# For each ASV in the given fly sample, count the number of reads and the percentage of the total reads in the fly sample that the ASV accounts for

allFlyASVTables <- data.frame()

for (k in 1:length(validFamilies)) {
  family <- validFamilies[[k]]
  familyMembers <- family@member_list

  # Put each FamilyMember object in its own variable
  currentFly <- NULL

  for (i in 1:length(familyMembers)) {
      currentFly <- c(currentFly,familyMembers[[i]])
    	}

  # Put fly samples in a list to loop over:
  currentFlyList <- list(currentFly)
  
  # Remove any NULL elements in the list (in case there is only 1 fly sample)
  currentFlyList <- plyr::compact(currentFlyList)

  # Loop over the fly sample list and populate the "fly perspective" data frame:
  #   allFlyASVTables - for each fly sample, holds the ASV, its percentage of the total reads in the sample, and whether it's shared with manure
  
  for (i in 1:length(currentFlyList[[1]])) {
  
    current_Fly <- currentFlyList[[1]][[i]]
    
    # build fly sample's identifier (FamGroup.1 or FamGroup.2)
    FlyDesignation <- current_Fly@sample_designation
    FlyIdentifier <- paste(family@family_group, FlyDesignation, sep = "/")
	
	Manure_ASV_counts = CountASVs("Manure-Manure", abdTable)

    # merge current fly sample's ASV counts with manure's
    FlyTable <- merge(current_Fly@ASV_counts,
                         Manure_ASV_counts,
                         by = "ASV", all = TRUE)
    
    # Replace column names with the samples' sample_designations for easier plotting downstream:
    
      currentColName <- names(FlyTable)[2]
      setnames(FlyTable, old = currentColName, new = familyMembers[[currentColName]]@sample_designation)
    
    # Collapse table to figure out what is shared from fly sample's perspective - remove any ASVs that have a count of NA in the current sample
    FlyTable <- FlyTable[!(is.na(FlyTable[[FlyDesignation]])), ]
    
    # Add a column to indicate if each ASV (row) in the collapsed table is shared with manure or not. If the ASV has a value of NA in manure column then it cannot be shared between current fly sample and manure
    
    FlyTable$shared <- ifelse(!(is.na(FlyTable$"Manure-Manure")), yes = TRUE, no = FALSE)
    FlyTable$sample <- FlyIdentifier

    # Prepare to add this table to a bigger list, but don't add until end!
    setnames(FlyTable, old = FlyDesignation, new = "read_count")
      
    # For this fly sample, calculate percent of reads that are from ASVs that are shared with manure, and not shared with manure (keep this table).
    FlySum <- sum(FlyTable$read_count)
    FlySharedSum <- sum(FlyTable$read_count[FlyTable$shared == TRUE])
    FlyNotSharedSum <- sum(FlyTable$read_count[FlyTable$shared == FALSE])
      
    # Add the flyTable (the one with ASVs as rows) to larger data frame
    FlyTablePercentages <- FlyTable %>%
      mutate("Fly_percent" = read_count/FlySum) %>%
      dplyr::select(-c("Manure-Manure"))
        
    allFlyASVTables <- rbind(allFlyASVTables, FlyTablePercentages)
  }
}

# Generate a plot that is similar to the one above, but by richness instead of relative abundance of reads.

# Group allFlyASVTables by sample and shared (T/F) then count

FlyPerspectiveRichness <- allFlyASVTables
array <- unlist(strsplit(FlyPerspectiveRichness$sample, "[/]"))
date <- array[c(TRUE, FALSE)]
FlyPerspectiveRichness$date <- date

FlyPerspectiveRichness.1 <- FlyPerspectiveRichness %>%
  group_by(date, sample, shared) %>%
  summarise(number_of_ASVs = n())
  
lst <- strsplit(FlyPerspectiveRichness.1$sample,'/')
v2 <- lapply(lst, `[`, 2)
v2 = substr(v2,1,4)
FlyPerspectiveRichness.1$sample.type <- v2

FlySharedOTUs <- FlyPerspectiveRichness.1
FlySharedOTUs <- filter(FlySharedOTUs, shared == TRUE)
FlySharedOTUs <- subset(FlySharedOTUs, select = c(date,sample,number_of_ASVs, sample.type))
FlySharedOTUs$sample.type <- factor(FlySharedOTUs$sample.type, levels=c("Endo", "Ecto"))
kruskal.test(number_of_ASVs~as.factor(sample.type),data=FlySharedOTUs) #Significant (p < 0.0001)
means <- aggregate(number_of_ASVs ~ sample.type, data = FlySharedOTUs, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(number_of_ASVs ~ sample.type, data = FlySharedOTUs, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="sample.type")
summary$sample.type <- factor(summary$sample.type, levels=c("Endo", "Ecto"))

ggplot(summary, aes(x = sample.type, y = number_of_ASVs.x)) +
    geom_errorbar(aes(ymin=number_of_ASVs.x-number_of_ASVs.y, ymax=number_of_ASVs.x+number_of_ASVs.y), width=.1) +
    geom_point()
  
#### Fig 5B (left, Arlington samples)

FlyPerspectiveRA <- allFlyASVTables
array <- unlist(strsplit(FlyPerspectiveRA$sample, "[/]"))
date <- array[c(TRUE, FALSE)]
FlyPerspectiveRA$date <- date

FlyPerspectiveRA.1 <- FlyPerspectiveRA %>%
  group_by(date, sample, shared) %>%
  summarise(tot_shared_reads = sum(read_count))
  
FlyPerspectiveRA.2 <- FlyPerspectiveRA.1 %>%
  group_by(sample) %>%
  summarise(tot_reads = sum(tot_shared_reads))

FlyPerspectiveRA.3 <- merge(x = FlyPerspectiveRA.1,y = FlyPerspectiveRA.2,by.x = "sample",by.y = "sample", all = T)

FlyPerspectiveRA.3$percent <- FlyPerspectiveRA.3$tot_shared_reads/FlyPerspectiveRA.3$tot_reads

#write.csv(F1PerspectiveRA.3, file = "ManurePerspectiveRAFlies.csv", quote = FALSE, row.names = FALSE)

data=read.csv("DCCManurePerspectiveRAFliesFormatted.csv")
kruskal.test(percent_shared~as.factor(fly.source),data=data) #Significant (p = 0.04489)
kruskal.test(percent_shared~as.factor(date),data=data) #Significant (p = 0.01791)
pairwise.wilcox.test(data$percent_shared, as.factor(data$date),
                 p.adjust.method = "BH", exact = FALSE) #Only significant: 7.9.2021 vs. 8.13.2021

data_summary <- aggregate(percent_shared ~ as.factor(fly.source), data,       # Create summary data
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("Endo", "Ecto"))
ggplot(data_summary) + #Fig. 5B
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, linewidth=1.3) +
    ylim(0,1)

#### Fig 5C (left, Arlington samples)

ps <- readRDS("ps_FieldWork2021_AJS_Final.rds")
sample_data(ps)$merge_factor <- c("Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","ArlM1-7.9.21","ArlM10-7.16.21","ArlM11-7.23.21","ArlM12-7.23.21","ArlM13-7.23.21","ArlM14-7.23.21","ArlM15-7.23.21","ArlM16-7.30.21","ArlM17-7.30.21","ArlM18-7.30.21","ArlM19-7.30.21","ArlM2-7.9.21","ArlM20-7.30.21","ArlM21-7.30.21","ArlM22-8.6.21","ArlM23-8.6.21","ArlM24-8.6.21","ArlM25-8.6.21","ArlM26-8.6.21","ArlM27-8.6.21","ArlM28-8.13.21","ArlM29-8.13.21","ArlM3-7.9.21","ArlM30-8.13.21","ArlM31-8.13.21","ArlM32-8.13.21","ArlM33-8.13.21","ArlM34-8.20.21","ArlM35-8.20.21","ArlM36-8.20.21","ArlM37-8.20.21","ArlM38-8.20.21","ArlM39-8.20.21","ArlM4-7.9.21","ArlM40-8.27.21","ArlM41-8.27.21","ArlM42-8.27.21","ArlM43-8.27.21","ArlM44-8.27.21","ArlM45-8.27.21","ArlM46-9.10.21","ArlM47-9.10.21","ArlM48-9.10.21","ArlM49-9.10.21","ArlM5-7.9.21","ArlM50-9.10.21","ArlM51-9.10.21","ArlM6-7.16.21","ArlM7-7.16.21","ArlM8-7.16.21","ArlM9-7.16.21","DCCM100-9.15.21","DCCM101-9.15.21","DCCM102-7.8.21","DCCM103-7.8.21","DCCM104-7.8.21","DCCM105-7.8.21","DCCM106-7.8.21","DCCM52-7.15.21","DCCM53-7.15.21","DCCM54-7.15.21","DCCM55-7.15.21","DCCM56-7.15.21","DCCM57-7.22.21","DCCM58-7.22.21","DCCM59-7.22.21","DCCM60-7.22.21","DCCM61-7.22.21","DCCM62-7.29.21","DCCM63-7.29.21","DCCM64-7.29.21","DCCM65-7.29.21","DCCM66-7.29.21","DCCM67-8.5.21","DCCM68-8.5.21","DCCM69-8.5.21","DCCM70-8.5.21","DCCM71-8.5.21","DCCM72-8.12.21","DCCM73-8.12.21","DCCM74-8.12.21","DCCM75-8.12.21","DCCM76-8.12.21","DCCM77-8.19.21","DCCM78-8.19.21","DCCM79-8.19.21","DCCM80-8.19.21","DCCM81-8.19.21","DCCM82-8.26.21","DCCM83-8.26.21","DCCM84-8.26.21","DCCM85-8.26.21","DCCM86-8.26.21","DCCM87-9.2.21","DCCM88-9.2.21","DCCM89-9.2.21","DCCM90-9.2.21","DCCM91-9.2.21","DCCM92-9.9.21","DCCM93-9.9.21","DCCM94-9.9.21","DCCM95-9.9.21","DCCM96-9.9.21","DCCM97-9.15.21","DCCM98-9.15.21","DCCM99-9.15.21","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly")
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

# Loop over the rows of the sample data data frame, building the objects by familiy group number.  But we only want to build it once, so when we run across the family group number again we don't want to re-build it.

for (i in 1:55) {
  # Get the current family group from the sample data table
  currentFamilyGroup <- sampleData[i, "sampling.date"]
  
  # Get all the current family names in familyClassList to check if 
  #   currentFamilyGroup has already been built:
  familyNames <- names(familyClassList)
  
  # If the current family group has not been built, build it:
  if (!(currentFamilyGroup %in% familyNames)){
    # Populate the family_group slot:
    # Create a list of sample names for all the samples in the current family group by subsetting the sample data df and extracting the sample names as a char vector.  Then, add to the memberList list that will be looped over to build each FamilyMember object for the current Family.
    members <- sampleData %>%
      filter(sampling.date == currentFamilyGroup) %>%
      dplyr::select(all_of("sample_id"))
    memberList <- as.list(members$sample_id)
    
    # Create an empty list to hold FamilyMember objects which will be added to the Family object in the end.
    currentFamilyMembersList <- list()
    
    # Loop over the sample names in memberList and build a new FamilyMember object for each sample:
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
    
    # Now that the member list is complete, get a total read sum to be put into the Family object (sum of each family member's total read count).
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

# A valid Family must have at least one manure sample and one fly sample:
# So a Family object should have a family_member list equal to length 3 or equal to length 2 where one of the sample_designation values is "Manure"

# Create a list to hold the valid Family objects:
validFamilies <- familyClassList

# Goal: generate a data frame holding all manure samples and all their ASVs (as rows). For each ASV in the given manure sample, count the number of reads and the percentage of the total reads in the manure sample that the ASV accounts for.

allManureASVTables <- data.frame()

for (k in 1:length(validFamilies)) {
  family <- validFamilies[[k]]
  familyMembers <- family@member_list

  # Put each FamilyMember object in its own variable
  currentManure <- NULL

  for (i in 1:length(familyMembers)) {
      currentManure <- c(currentManure,familyMembers[[i]])
    	}

  # Put manure samples in a list to loop over:
  currentManureList <- list(currentManure)
  
  # Remove any NULL elements in the list (in case there is only 1 fly sample)
  currentManureList <- plyr::compact(currentManureList)

  # Loop over the manure sample list and populate the "manure perspective" data frame:
  #   allManureASVTables - for each manure sample, holds the ASV, its percentage of the total reads in the manure sample, and whether it's shared with a fly sample
  
  for (i in 1:length(currentManureList[[1]])) {
  
    current_Manure <- currentManureList[[1]][[i]]
    
    # build manure sample's identifier (FamGroup.1 or FamGroup.2)
    ManureDesignation <- current_Manure@sample_designation
    ManureIdentifier <- paste(family@family_group, ManureDesignation, sep = "/")
	
	Fly_ASV_counts = CountASVs("Fly-Fly", abdTable)

    # merge current manure sample's ASV counts with flies'
    ManureTable <- merge(current_Manure@ASV_counts,
                         Fly_ASV_counts,
                         by = "ASV", all = TRUE)
    
    # Replace column names with the samples' sample_designations for easier plotting downstream:
    
      currentColName <- names(ManureTable)[2]
      setnames(ManureTable, old = currentColName, new = familyMembers[[currentColName]]@sample_designation)
    
    # Collapse table to figure out what is shared from manure sample's perspective - remove any ASVs that have a count of NA in the current manure sample
    ManureTable <- ManureTable[!(is.na(ManureTable[[ManureDesignation]])), ]
    
    # Add a column to indicate if each ASV (row) in the collapsed table is shared with a fly sample or not. If the ASV has a value of NA in fly column then it cannot be shared between current manure sample and flies.
    
    ManureTable$shared <- ifelse(!(is.na(ManureTable$"Fly-Fly")), yes = TRUE, no = FALSE)
    ManureTable$sample <- ManureIdentifier

    # Prepare to add this table to a bigger list, but don't add until end!
    setnames(ManureTable, old = ManureDesignation, new = "read_count")
      
    # For this manure sample, calculate percent of reads that are from ASVs that are shared with flies, and not shared with flies (keep this table).
    ManureSum <- sum(ManureTable$read_count)
    ManureSharedSum <- sum(ManureTable$read_count[ManureTable$shared == TRUE])
    ManureNotSharedSum <- sum(ManureTable$read_count[ManureTable$shared == FALSE])
      
    # Add the ManureTable (the one with ASVs as rows) to larger data frame
    ManureTablePercentages <- ManureTable %>%
      mutate("Manure_percent" = read_count/ManureSum) %>%
      dplyr::select(-c("Fly-Fly"))
    
    allManureASVTables <- rbind(allManureASVTables, ManureTablePercentages)
  }
}

ManurePerspectiveRA <- allManureASVTables
array <- unlist(strsplit(ManurePerspectiveRA$sample, "[/]"))
date <- array[c(TRUE, FALSE)]
ManurePerspectiveRA$date <- date

ManurePerspectiveRA.1 <- ManurePerspectiveRA %>%
  group_by(date, sample, shared) %>%
  summarise(tot_shared_reads = sum(read_count))
  
ManurePerspectiveRA.2 <- ManurePerspectiveRA.1 %>%
  group_by(sample) %>%
  summarise(tot_reads = sum(tot_shared_reads))

ManurePerspectiveRA.3 <- merge(x = ManurePerspectiveRA.1,y = ManurePerspectiveRA.2,by.x = "sample",by.y = "sample", all = T)

ManurePerspectiveRA.3$percent <- ManurePerspectiveRA.3$tot_shared_reads/ManurePerspectiveRA.3$tot_reads

#write.csv(ManurePerspectiveRA.3, file = "FlyPerspectiveRAManure.csv", quote = FALSE, row.names = FALSE)

data=read.csv("DCCFlyPerspectiveRAManureFormatted.csv")
kruskal.test(percent_shared~as.factor(date),data=data) #NS (p = 0.9071)

mean <- mean(data$percent_shared)
ci <- 1.96*(sd(data$percent_shared) / sqrt(length(data$percent_shared)))                              
ggplot(data=NULL) + #Fig. 5C
    geom_bar( aes(x="", y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x="", ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)
