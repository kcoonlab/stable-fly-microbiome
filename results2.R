set.seed(123456)

source("./Helper_Functions.R") #From CITATION
source("./TwinMom16S_Shared_Functions.R") #From CITATION

# Load required packages
library(phyloseq)
library(dplyr)
library(tidyverse)
library(data.table)
library(ggplot2)

#### Note: The below code is specific to Arlington-derived samples ####

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

#### Goal: Calculate the number of shared ASVs between internal/external fly samples and manure samples collected from the same facility (Fig. 5A) ####

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

#-------------------------------------------------------#

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

#-------------------------------------------------------#

CountASVs <- function(sampleName, abundanceTable) {
  sampleAbundance <- abundanceTable %>%
    dplyr::select(sampleName) %>%
    rownames_to_column(var = "ASV")
  sampleAbundanceDF <- filter(sampleAbundance, sampleAbundance[, sampleName] > 0)
  return(sampleAbundanceDF)
}

#-------------------------------------------------------#

GetSampleSum <- function(sampleName, physeqObject) {
  named_sum <- sample_sums(physeqObject)[sampleName]
  numeric_sum <- unname(named_sum)
}

#-------------------------------------------------------#

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
      dplyr::select(as.vector("sample_id"))
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
# A valid Family must have at least one manure and one fly sample:
# So a Family object should have a family_member list equal to length 3 or
# equal to length 2 where one of the sample_designation values is "Manure"

# Create a list to hold the valid Family objects:
validFamilies <- familyClassList

# Generate a data frame holding all fly samples and all their ASVs (as rows). For each
# ASV in the given fly sample, count the number of reads and
# the percentage of the total reads in the sample that the ASV accounts for

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

  # Loop over the fly sample list and populate the data frame:
  #   allFlyASVTables - for each fly sample, holds the ASV, 
  #     its percentage of the total reads in the sample, and whether it's shared with manure collected from the same facility
  
  for (i in 1:length(currentFlyList[[1]])) {
  
    current_Fly <- currentFlyList[[1]][[i]]
    #currentFlyName <- current_fly@name
    
    # build fly sample's identifier (FamGroup.1 or FamGroup.2)
    FlyDesignation <- current_Fly@sample_designation
    FlyIdentifier <- paste(family@family_group, FlyDesignation, sep = "/")
	
	Manure_ASV_counts = CountASVs("Manure-Manure", abdTable)

    # merge current fly sample's ASV counts with manure's
    FlyTable <- merge(current_Fly@ASV_counts,
                         Manure_ASV_counts,
                         by = "ASV", all = TRUE)
    
    # Replace column names with the samples' sample_designations for easier
    #   plotting downstream:
    
      currentColName <- names(FlyTable)[2]
      setnames(FlyTable, old = currentColName, new = familyMembers[[currentColName]]@sample_designation)
    
    # Collapse table to figure out what is shared from fly sample's perspective -
    #   remove any ASVs that have a count of NA in the current sample
    FlyTable <- FlyTable[!(is.na(FlyTable[[FlyDesignation]])), ]
    
    # Add a column to indicate if each ASV (row) in the collapsed table is 
    #   shared with manure or not. If the ASV has a value of NA in manure column
    #   then it cannot be shared between current fly sample and manure
    
    FlyTable$shared <- ifelse(!(is.na(FlyTable$"Manure-Manure")), yes = TRUE, no = FALSE)
    FlyTable$sample <- FlyIdentifier

    # Prepare to add this table to a bigger list, but don't add until end!
    setnames(FlyTable, old = FlyDesignation, new = "read_count")
      
    # For this fly sample, calculate percent of reads that are from ASVs that
    #   are shared with mom, and not shared with mom (keep this table).
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

#### Goal: Calculate relative abundance of fly-associated taxa that are shared with manure collected from the same facility (or not) (Fig. 5B) ####

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
FlyPerspectiveRA.3$source <- ifelse(str_detect(FlyPerspectiveRA.3$sample,"Endo"), "Endo", "Ecto")
FlyPerspectiveRA.3 <- FlyPerspectiveRA.3[FlyPerspectiveRA.3$shared == "TRUE",]		 
		
kruskal.test(percent~as.factor(source),data=FlyPerspectiveRA.3) #P = 0.04602
kruskal.test(percent~as.factor(date),data=FlyPerspectiveRA.3) #P = 0.01786
pairwise.wilcox.test(FlyPerspectiveRA.3$percent, as.factor(FlyPerspectiveRA.3$date),
                 p.adjust.method = "BH", exact = FALSE) #Only significant: 7.9.2021 vs. 8.13.2021 (P = 0.030)
		 
data_summary <- aggregate(percent ~ as.factor(source), FlyPerspectiveRA.3,
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("Endo", "Ecto"))
ggplot(data_summary) +
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)
		 	 
#### Goal: Calculate the number of shared ASVs between internal/external fly samples and manure samples collected from the same facility across different sampling dates (Fig. S7) ####
			  
FlySharedOTUs$date <- factor(FlySharedOTUs$date, levels=c("7.9.2021","7.16.2021","7.23.2021","7.30.2021","8.6.2021","8.13.2021","8.20.2021","8.27.2021","9.10.2021"))
FlySharedOTUs.Endo <- FlySharedOTUs[FlySharedOTUs$sample.type=="Endo",]
FlySharedOTUs.Ecto <- FlySharedOTUs[FlySharedOTUs$sample.type=="Ecto",]
kruskal.test(number_of_ASVs~as.factor(date),data=FlySharedOTUs.Endo) #NS (P = 0.5275)
kruskal.test(number_of_ASVs~as.factor(date),data=FlySharedOTUs.Ecto) #NS (P = 0.0781)
			  
means <- aggregate(number_of_ASVs ~ date, data = FlySharedOTUs.Ecto, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(number_of_ASVs ~ date, data = FlySharedOTUs.Ecto, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="date")
summary$date <- factor(summary$date, levels=c("7.9.2021","7.16.2021","7.23.2021","7.30.2021","8.6.2021","8.13.2021","8.20.2021","8.27.2021","9.10.2021"))
ggplot(summary, aes(x = date, y = number_of_ASVs.x)) +
    geom_errorbar(aes(ymin=number_of_ASVs.x-number_of_ASVs.y, ymax=number_of_ASVs.x+number_of_ASVs.y), width=.1) +
    geom_point()

#### Goal: Calculate the number of shared ASVs between internal/external fly samples and manure samples collected from the same facility across different sampling locations (Fig. S8) ####

means <- aggregate(number_of_ASVs ~ trap, data = FlySharedOTUs.Endo, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(number_of_ASVs ~ trap, data = FlySharedOTUs.Endo, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="trap")
summary$date <- factor(summary$trap, levels=c("B1-East", "B1-West", "B2-East", "B2-West"))
ggplot(summary, aes(x = trap, y = number_of_ASVs.x)) +
    geom_errorbar(aes(ymin=number_of_ASVs.x-number_of_ASVs.y, ymax=number_of_ASVs.x+number_of_ASVs.y), width=.1) +
    geom_point()

means <- aggregate(number_of_ASVs ~ trap, data = FlySharedOTUs.Ecto, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(number_of_ASVs ~ trap, data = FlySharedOTUs.Ecto, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="trap")
summary$date <- factor(summary$trap, levels=c("B1-East", "B1-West", "B2-East", "B2-West"))
ggplot(summary, aes(x = date, y = number_of_ASVs.x)) +
    geom_errorbar(aes(ymin=number_of_ASVs.x-number_of_ASVs.y, ymax=number_of_ASVs.x+number_of_ASVs.y), width=.1) +
    geom_point()

#### Goal: Calculate relative abundance of fly-associated taxa that are shared with manure collected from the same facility (or not), across different sampling dates (Fig. S9) ####

data.endo <- data[data$fly.source=="Endo",]
data_summary <- aggregate(percent_shared ~ as.factor(date), data.endo,
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("7.9.2021","7.16.2021","7.23.2021","7.30.2021","8.6.2021","8.13.2021","8.20.2021","8.27.2021","9.10.2021"))
ggplot(data_summary) + #Fig. S9
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)

data.ecto <- data[data$fly.source=="Ecto",]
data_summary <- aggregate(percent_shared ~ as.factor(date), data.ecto,
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("7.9.2021","7.16.2021","7.23.2021","7.30.2021","8.6.2021","8.13.2021","8.20.2021","8.27.2021","9.10.2021"))
ggplot(data_summary) + #Fig. S9
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)

#### Goal: Identify ASVs that are commonly shared between fly and manure samples collected from the same facility (Fig. 6) ####
			  
ps.Merged.manure <- subset_taxa(ps.Merged,otu_table(ps.Merged)["Manure-Manure",]>0)
ps.Merged.manure <- subset_samples(ps.Merged.manure, pooltype != "Manure")
asvfac <- rownames(otu_table(ps.Merged.manure))
asvtab = apply(otu_table(ps.Merged.manure), MARGIN = 2, function(x) {
    tapply(x, INDEX = asvfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
asvtab <- t(asvtab)
asvtab <- rowSums(asvtab > 0)
df.asvtab <- data.frame(asvtab)
df.asvtab.out <- df.asvtab %>%
  filter(asvtab > 34) %>%
  arrange(asvtab)
df.asvtab.out$asv <- rownames(df.asvtab.out)
df.asvtab.out                   
tax_table(ps.Merged.manure)[df.asvtab.out$asv,]                        

ggplot(df.asvtab.out, aes(x = reorder(asv, -asvtab,sum), y = asvtab)) +
  geom_bar(fill = "#DCDCDC", stat = "identity")  

#### Goal: Plot relative abundance of commonly shared ASVs in internal/external fly vs. manure samples collected from the same facility (Fig. 7) ####
                          
ps.Merged <- merge_samples(ps.arlington,"sample.type")
ps.Merged <- transform_sample_counts(ps.Merged, function(x) x / sum(x) )
row_names <- row.names(sample_data(ps.Merged))
sample_data(ps.Merged)$sample_id <- row_names
otus <- otu_table(ps.Merged)
otus <- as.data.frame(t(otus))
row_names <- row.names(otus)
otus$ASV <- row_names
                                     
otus <- otus[otus$ASV %in% c("d56a560e43e7ab41839a1c8a10f65d36", "aaa9219facfcb711b04e0d3c190b571a", "3289d95c4a81c516faa9d263f55202d2", "8da7704aad448995f55464ee650a4452", "66402d8adb61bb82bd8d7fb58d832c10",
                "6207db2806b88c7f75541b8de1c42f50", "91f2efe0637f2f2062c1252c32a8e4eb", "350bce411d627b38ef18aeaebd7d4f69", "dcfd8312bf9f39a54aab3707c3b140df", "945bb343f8dc3dd656b9ca15c3c0f09d",
                "e0597b0224419e26dd82616c75f8c1d2", "bcaaa3d15ff438cfa14dbec5d86b0eb8", "5dbd59fa8200705770f4d068d6572dc3", "b77e07971c99cc0a36c99bd7cad83c16", "b40dfe8590cdabab904a9931210d2805",
                "ff27974c65664de75e6421a62a0448a3", "15d888a8d8e6c158974710b615dcbcdf"), ]

otus$ASV <- factor(otus$ASV,levels = c("d56a560e43e7ab41839a1c8a10f65d36", "aaa9219facfcb711b04e0d3c190b571a", "3289d95c4a81c516faa9d263f55202d2", "8da7704aad448995f55464ee650a4452", "66402d8adb61bb82bd8d7fb58d832c10",
                "6207db2806b88c7f75541b8de1c42f50", "91f2efe0637f2f2062c1252c32a8e4eb", "350bce411d627b38ef18aeaebd7d4f69", "dcfd8312bf9f39a54aab3707c3b140df", "945bb343f8dc3dd656b9ca15c3c0f09d",
                "e0597b0224419e26dd82616c75f8c1d2", "bcaaa3d15ff438cfa14dbec5d86b0eb8", "5dbd59fa8200705770f4d068d6572dc3", "b77e07971c99cc0a36c99bd7cad83c16", "b40dfe8590cdabab904a9931210d2805",
                "ff27974c65664de75e6421a62a0448a3", "15d888a8d8e6c158974710b615dcbcdf"))

ASVNames <- c("d56a560e43e7ab41839a1c8a10f65d36", "aaa9219facfcb711b04e0d3c190b571a", "3289d95c4a81c516faa9d263f55202d2", "8da7704aad448995f55464ee650a4452", "66402d8adb61bb82bd8d7fb58d832c10",
                "6207db2806b88c7f75541b8de1c42f50", "91f2efe0637f2f2062c1252c32a8e4eb", "350bce411d627b38ef18aeaebd7d4f69", "dcfd8312bf9f39a54aab3707c3b140df", "945bb343f8dc3dd656b9ca15c3c0f09d",
                "e0597b0224419e26dd82616c75f8c1d2", "bcaaa3d15ff438cfa14dbec5d86b0eb8", "5dbd59fa8200705770f4d068d6572dc3", "b77e07971c99cc0a36c99bd7cad83c16", "b40dfe8590cdabab904a9931210d2805",
                "ff27974c65664de75e6421a62a0448a3", "15d888a8d8e6c158974710b615dcbcdf")

myColors <- c("#edf3f7", "#d2e2ef", "#9fc1dc", "#4292c6", "#2171b5", "#fee6ce", "#fdae6b", "#fd8d3c", "#f16913", "#d94801", "#a63603", 
              "#7f2704", "#f1efe1", "#e9e9e7", "#d4d4d4", "#bababa", "#878787")
names(myColors) <- ASVNames

#Manure
ggplot(otus, aes(fill=ASV, y=Manure, x="")) + #Fig. 7
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)
				     
#Ecto
ggplot(otus, aes(fill=ASV, y=Ecto, x="")) + #Fig. 7
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)
				     
#Endo
ggplot(otus, aes(fill=ASV, y=Ecto, x="")) + #Fig. 7
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)

#### Goal: Plot relative abundance of commonly shared ASVs in internal/external fly vs. manure samples collected from the same facility, across different sampling dates (Fig. S11) ####

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
                        
#-------------------------------------------------------#

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

#-------------------------------------------------------#

CountASVs <- function(sampleName, abundanceTable) {
  sampleAbundance <- abundanceTable %>%
    dplyr::select(sampleName) %>%
    rownames_to_column(var = "ASV")
  sampleAbundanceDF <- filter(sampleAbundance, sampleAbundance[, sampleName] > 0)
  return(sampleAbundanceDF)
}

#-------------------------------------------------------#

GetSampleSum <- function(sampleName, physeqObject) {
  named_sum <- sample_sums(physeqObject)[sampleName]
  numeric_sum <- unname(named_sum)
}

#-------------------------------------------------------#

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
      dplyr::select(as.vector("sample_id"))
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
  currentManure <- NULL
  currentEndo <- NULL
  currentEcto <- NULL

  for (i in 1:length(familyMembers)) {
    if (str_starts(familyMembers[[i]]@sample_designation, "Manure")) {
      currentManure <- familyMembers[[i]]
    } else if (str_starts(familyMembers[[i]]@sample_designation, "Endo")) {
      currentEndo <- familyMembers[[i]]
    } else if (str_starts(familyMembers[[i]]@sample_designation, "Ecto")) {
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

allF1ASVTables <- allF1ASVTables[allF1ASVTables$ASV %in% c("d56a560e43e7ab41839a1c8a10f65d36", "aaa9219facfcb711b04e0d3c190b571a", "3289d95c4a81c516faa9d263f55202d2", "8da7704aad448995f55464ee650a4452", "66402d8adb61bb82bd8d7fb58d832c10",
                "6207db2806b88c7f75541b8de1c42f50", "91f2efe0637f2f2062c1252c32a8e4eb", "350bce411d627b38ef18aeaebd7d4f69", "dcfd8312bf9f39a54aab3707c3b140df", "945bb343f8dc3dd656b9ca15c3c0f09d",
                "e0597b0224419e26dd82616c75f8c1d2", "bcaaa3d15ff438cfa14dbec5d86b0eb8", "5dbd59fa8200705770f4d068d6572dc3", "b77e07971c99cc0a36c99bd7cad83c16", "b40dfe8590cdabab904a9931210d2805",
                "ff27974c65664de75e6421a62a0448a3", "15d888a8d8e6c158974710b615dcbcdf"), ]

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
ASVTax <- ASVTax[ASVTax$ASV %in% c("d56a560e43e7ab41839a1c8a10f65d36", "aaa9219facfcb711b04e0d3c190b571a", "3289d95c4a81c516faa9d263f55202d2", "8da7704aad448995f55464ee650a4452", "66402d8adb61bb82bd8d7fb58d832c10",
                "6207db2806b88c7f75541b8de1c42f50", "91f2efe0637f2f2062c1252c32a8e4eb", "350bce411d627b38ef18aeaebd7d4f69", "dcfd8312bf9f39a54aab3707c3b140df", "945bb343f8dc3dd656b9ca15c3c0f09d",
                "e0597b0224419e26dd82616c75f8c1d2", "bcaaa3d15ff438cfa14dbec5d86b0eb8", "5dbd59fa8200705770f4d068d6572dc3", "b77e07971c99cc0a36c99bd7cad83c16", "b40dfe8590cdabab904a9931210d2805",
                "ff27974c65664de75e6421a62a0448a3", "15d888a8d8e6c158974710b615dcbcdf"), ]
                          
allF1ASVTables$ASV <- factor(allF1ASVTables$ASV,levels = c("d56a560e43e7ab41839a1c8a10f65d36", "aaa9219facfcb711b04e0d3c190b571a", "3289d95c4a81c516faa9d263f55202d2", "8da7704aad448995f55464ee650a4452", "66402d8adb61bb82bd8d7fb58d832c10",
                "6207db2806b88c7f75541b8de1c42f50", "91f2efe0637f2f2062c1252c32a8e4eb", "350bce411d627b38ef18aeaebd7d4f69", "dcfd8312bf9f39a54aab3707c3b140df", "945bb343f8dc3dd656b9ca15c3c0f09d",
                "e0597b0224419e26dd82616c75f8c1d2", "bcaaa3d15ff438cfa14dbec5d86b0eb8", "5dbd59fa8200705770f4d068d6572dc3", "b77e07971c99cc0a36c99bd7cad83c16", "b40dfe8590cdabab904a9931210d2805",
                "ff27974c65664de75e6421a62a0448a3", "15d888a8d8e6c158974710b615dcbcdf"))

ASVNames <- c("d56a560e43e7ab41839a1c8a10f65d36", "aaa9219facfcb711b04e0d3c190b571a", "3289d95c4a81c516faa9d263f55202d2", "8da7704aad448995f55464ee650a4452", "66402d8adb61bb82bd8d7fb58d832c10",
                "6207db2806b88c7f75541b8de1c42f50", "91f2efe0637f2f2062c1252c32a8e4eb", "350bce411d627b38ef18aeaebd7d4f69", "dcfd8312bf9f39a54aab3707c3b140df", "945bb343f8dc3dd656b9ca15c3c0f09d",
                "e0597b0224419e26dd82616c75f8c1d2", "bcaaa3d15ff438cfa14dbec5d86b0eb8", "5dbd59fa8200705770f4d068d6572dc3", "b77e07971c99cc0a36c99bd7cad83c16", "b40dfe8590cdabab904a9931210d2805",
                "ff27974c65664de75e6421a62a0448a3", "15d888a8d8e6c158974710b615dcbcdf")

myColors <- c("#edf3f7", "#d2e2ef", "#9fc1dc", "#4292c6", "#2171b5", "#fee6ce", "#fdae6b", "#fd8d3c", "#f16913", "#d94801",
        "#a63603", "#7f2704", "#f1efe1", "#e9e9e7", "#d4d4d4", "#bababa", "#878787")
names(myColors) <- ASVNames

F1PerspectivePlot <- ggplot(allF1ASVTables,
                                aes(x = sample, y = F1_percent, 
                                    fill = ASV)) +
  geom_bar(stat = "identity") +
scale_fill_manual(name = "ASV", values = myColors)
F1PerspectivePlot

#### Goal: Calculate relative abundance of manure-associated taxa that are shared with flies collected from the same facility (or not) (Fig. 5C) ####

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

#-------------------------------------------------------#

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

#-------------------------------------------------------#

CountASVs <- function(sampleName, abundanceTable) {
  sampleAbundance <- abundanceTable %>%
    dplyr::select(sampleName) %>%
    rownames_to_column(var = "ASV")
  sampleAbundanceDF <- filter(sampleAbundance, sampleAbundance[, sampleName] > 0)
  return(sampleAbundanceDF)
}

#-------------------------------------------------------#

GetSampleSum <- function(sampleName, physeqObject) {
  named_sum <- sample_sums(physeqObject)[sampleName]
  numeric_sum <- unname(named_sum)
}

#-------------------------------------------------------#

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

for (i in 1:51) {
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
      dplyr::select(as.vector("sample_id"))
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
  currentManure <- NULL

  for (i in 1:length(familyMembers)) {
      currentManure <- c(currentManure,familyMembers[[i]])
    	}

  # Put infants in a list to loop over:
  currentF1List <- list(currentManure)
  
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
	
	Fly_ASV_counts = CountASVs("Fly-Fly", abdTable)

    # merge current infant's ASV counts with mother's
    F1Table <- merge(current_F1@ASV_counts,
                         Fly_ASV_counts,
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
    
    F1Table$shared <- ifelse(!(is.na(F1Table$"Fly-Fly")), yes = TRUE, no = FALSE)
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
      dplyr::select(-c("Fly-Fly"))
    
    # Add infant sample name to table
    #infantTablePercentages <- cbind(infantTablePercentages,
                                    #"infant_sample" = currentInfantName)
    
    allF1ASVTables <- rbind(allF1ASVTables, F1TablePercentages)
  }
}

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
F1PerspectiveRA.3 <- F1PerspectiveRA.3[F1PerspectiveRA.3$shared == "TRUE",]	
				     
kruskal.test(percent~as.factor(date),data=F1PerspectiveRA.3) #NS (p = 1)

mean <- mean(F1PerspectiveRA.3$percent)
ci <- 1.96*(sd(F1PerspectiveRA.3$percent) / sqrt(length(F1PerspectiveRA.3$percent)))        

ggplot(data=NULL) +
    geom_bar( aes(x="", y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x="", ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)

#### Goal: Calculate relative abundance of manure-associated taxa that are shared with flies collected from the same facility (or not), across different sampling dates (Fig. S10) ####
				     
means <- aggregate(percent ~ date, data = F1PerspectiveRA.3, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(percent ~ date, data = F1PerspectiveRA.3, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="date")
summary$date <- factor(summary$date, levels=c("7.9.21","7.16.21","7.23.21","7.30.21","8.6.21","8.13.21","8.20.21","8.27.21","9.10.21"))

ggplot(summary, aes(x = date, y = percent.x)) +
    geom_errorbar(aes(ymin=percent.x-percent.y, ymax=percent.x+percent.y), width=.1) +
    geom_point()

#### Note: The below code is specific to DCC-derived samples ####

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

#### Goal: Calculate the number of shared ASVs between internal/external fly samples and manure samples collected from the same facility (Fig. 5A) ####

sample_data(ps)$merge_factor <- c("Endo10A-7.16.2021","Endo10B-7.16.2021","Endo10D-7.16.2021","Endo11A-7.22.2021","Endo11B-7.22.2021","Endo11D-7.22.2021","Endo11E-7.22.2021","Endo12B-7.23.2021","Endo12C-7.23.2021","Endo12D-7.23.2021","Endo12E-7.23.2021","Endo13A-7.23.2021","Endo13B-7.23.2021","Endo13C-7.23.2021","Endo13D-7.23.2021","Endo13E-7.23.2021","Endo14A-7.23.2021","Endo14B-7.23.2021","Endo14C-7.23.2021","Endo15A-7.23.2021","Endo15B-7.23.2021","Endo15C-7.23.2021","Endo15D-7.23.2021","Endo16A-7.29.2021","Endo17C-7.30.2021","Endo18A-7.30.2021","Endo18B-7.30.2021","Endo19A-7.30.2021","Endo19B-7.30.2021","Endo1A-7.8.2021","Endo1B-7.8.2021","Endo20A-7.30.2021","Endo20B-7.30.2021","Endo21A-8.5.2021","Endo22A-8.6.2021","Endo22B-8.6.2021","Endo22C-8.6.2021","Endo22D-8.6.2021","Endo22E-8.6.2021","Endo23A-8.6.2021","Endo23B-8.6.2021","Endo23C-8.6.2021","Endo23D-8.6.2021","Endo23E-8.6.2021","Endo24A-8.6.2021","Endo24B-8.6.2021","Endo24C-8.6.2021","Endo24D-8.6.2021","Endo24E-8.6.2021","Endo25A-8.6.2021","Endo25B-8.6.2021","Endo25C-8.6.2021","Endo25D-8.6.2021","Endo27A-8.13.2021","Endo27B-8.13.2021","Endo27C-8.13.2021","Endo27D-8.13.2021","Endo27E-8.13.2021","Endo28A-8.13.2021","Endo28B-8.13.2021","Endo29A-8.13.2021","Endo29B-8.13.2021","Endo2A-7.9.2021","Endo2B-7.9.2021","Endo30A-8.13.2021","Endo30B-8.13.2021","Endo30C-8.13.2021","Endo30D-8.13.2021","Endo32A-8.20.2021","Endo32D-8.20.2021","Endo33A-8.20.2021","Endo33B-8.20.2021","Endo34A-8.20.2021","Endo34B-8.20.2021","Endo34C-8.20.2021","Endo35A-8.20.2021","Endo37A-8.27.2021","Endo37B-8.27.2021","Endo37C-8.27.2021","Endo37D-8.27.2021","Endo37E-8.27.2021","Endo38A-8.27.2021","Endo39A-8.27.2021","Endo39B-8.27.2021","Endo39C-8.27.2021","Endo39D-8.27.2021","Endo3A-7.9.2021","Endo3C-7.9.2021","Endo3D-7.9.2021","Endo3E-7.9.2021","Endo40A-8.27.2021","Endo40B-8.27.2021","Endo40C-8.27.2021","Endo41A-9.9.2021","Endo42A-9.10.2021","Endo42B-9.10.2021","Endo42C-9.10.2021","Endo42D-9.10.2021","Endo42E-9.10.2021","Endo44A-9.10.2021","Endo44B-9.10.2021","Endo44C-9.10.2021","Endo44D-9.10.2021","Endo44E-9.10.2021","Endo46A-9.15.2021","Endo4A-7.9.2021","Endo4B-7.9.2021","Endo4C-7.9.2021","Endo4D-7.9.2021","Endo5A-7.9.2021","Endo5B-7.9.2021","Endo6A-7.15.2021","Endo6B-7.15.2021","Endo6C-7.15.2021","Endo6E-7.15.2021","Endo7B-7.16.2021","Endo7C-7.16.2021","Endo7D-7.16.2021","Endo7E-7.16.2021","Endo8A-7.16.2021","Endo8C-7.16.2021","Endo8D-7.16.2021","Endo8E-7.16.2021","Endo9A-7.16.2021","Endo9B-7.16.2021","Endo9C-7.16.2021","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Manure-Manure","Ecto24-8.6.2021","Ecto10-7.16.2021","Ecto11-7.22.2021","Ecto12-7.23.2021","Ecto13-7.23.2021","Ecto14-7.23.2021","Ecto16-7.29.2021","Ecto17-7.30.2021","Ecto18-7.30.2021","Ecto19-7.30.2021","Ecto20-7.30.2021","Ecto21-8.5.2021","Ecto22-8.6.2021","Ecto23-8.6.2021","Ecto25-8.6.2021","Ecto27-8.13.2021","Ecto28-8.13.2021","Ecto29-8.13.2021","Ecto30-8.13.2021","Ecto32-8.20.2021","Ecto33-8.20.2021","Ecto34-8.20.2021","Ecto35-8.20.2021","Ecto37-8.27.2021","Ecto38-8.27.2021","Ecto39-8.27.2021","Ecto40-8.27.2021","Ecto41-9.9.2021","Ecto42-9.10.2021","Ecto44-9.10.2021","Ecto7-7.16.2021","Ecto9-7.16.2021")
ps.arlington <- subset_samples(ps, location=="DCC")
ps.Merged <- merge_samples(ps.arlington,"merge_factor")
row_names <- row.names(sample_data(ps.Merged))
sample_data(ps.Merged)$sample_id <- row_names
lst <- strsplit(row_names,'-')
v1 <- lapply(lst, `[`, 1)
v2 <- lapply(lst, `[`, 2)
sample_data(ps.Merged)$pooltype <- v1
sample_data(ps.Merged)$sampling.date <- v2

#-------------------------------------------------------#

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

#-------------------------------------------------------#

CountASVs <- function(sampleName, abundanceTable) {
  sampleAbundance <- abundanceTable %>%
    dplyr::select(sampleName) %>%
    rownames_to_column(var = "ASV")
  sampleAbundanceDF <- filter(sampleAbundance, sampleAbundance[, sampleName] > 0)
  return(sampleAbundanceDF)
}

#-------------------------------------------------------#

GetSampleSum <- function(sampleName, physeqObject) {
  named_sum <- sample_sums(physeqObject)[sampleName]
  numeric_sum <- unname(named_sum)
}

#-------------------------------------------------------#

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
      dplyr::select(as.vector("sample_id"))
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
kruskal.test(number_of_ASVs~as.factor(sample.type),data=F1SharedOTUs) #NS (p = 0.6703)
means <- aggregate(number_of_ASVs ~ sample.type, data = F1SharedOTUs, 
          FUN = function(x) c(mean = mean(x)))

cis <- aggregate(number_of_ASVs ~ sample.type, data = F1SharedOTUs, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))

summary <- merge(means,cis,by="sample.type")
summary$sample.type <- factor(summary$sample.type, levels=c("Endo", "Ecto"))

ggplot(summary, aes(x = sample.type, y = number_of_ASVs.x)) + #Fig. 5A
    geom_errorbar(aes(ymin=number_of_ASVs.x-number_of_ASVs.y, ymax=number_of_ASVs.x+number_of_ASVs.y), width=.1) +
    geom_point()
		 
#### Goal: Calculate relative abundance of fly-associated taxa that are shared with manure collected from the same facility (or not) (Fig. 5B) ####

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
            file = "DCCManurePerspectiveRAFlies.csv",
            quote = FALSE, row.names = FALSE)

data=read.csv("DCCManurePerspectiveRAFliesFormatted.csv")
kruskal.test(percent_shared~as.factor(fly.source),data=data) #NS (p = 0.8318)
kruskal.test(percent_shared~as.factor(date),data=data) #NS (p = 0.7612)

data_summary <- aggregate(percent_shared ~ as.factor(fly.source), data,       # Create summary data
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("Endo", "Ecto"))
ggplot(data_summary) + #Fig. 5B
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)
			  
#### Goal: Calculate the number of shared ASVs between internal/external fly samples and manure samples collected from the same facility across different sampling dates (Fig. S7) ####

F1SharedOTUs$date <- factor(F1SharedOTUs$date, levels=c("7.8.2021","7.15.2021","7.22.2021","7.29.2021","8.5.2021","9.9.2021","9.15.2021"))
F1SharedOTUs.Endo <- F1SharedOTUs[F1SharedOTUs$sample.type=="Endo",]
F1SharedOTUs.Ecto <- F1SharedOTUs[F1SharedOTUs$sample.type=="Ecto",]
kruskal.test(number_of_ASVs~as.factor(date),data=F1SharedOTUs.Endo) #NS (p = 0.2037)
means <- aggregate(number_of_ASVs ~ date, data = F1SharedOTUs.Endo, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(number_of_ASVs ~ date, data = F1SharedOTUs.Endo, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="date")
summary$date <- factor(summary$date, levels=c("7.8.2021","7.15.2021","7.22.2021","7.29.2021","8.5.2021","9.9.2021","9.15.2021"))
ggplot(summary, aes(x = date, y = number_of_ASVs.x)) + #Fig. S7 (right)
    geom_errorbar(aes(ymin=number_of_ASVs.x-number_of_ASVs.y, ymax=number_of_ASVs.x+number_of_ASVs.y), width=.1) + 
    geom_point()
ggplot(F1SharedOTUs.Ecto, aes(x=date, y=number_of_ASVs)) + geom_point() # Fig. S7 (right)
		 
#### Goal: Calculate the number of shared ASVs between internal/external fly samples and manure samples collected from the same facility across different sampling locations (Fig. S8) ####
#### Goal: Calculate relative abundance of fly-associated taxa that are shared with manure collected from the same facility (or not), across different sampling dates (Fig. S9) ####

data.endo <- data[data$fly.source=="Endo",]
data_summary <- aggregate(percent_shared ~ as.factor(date), data.endo,       # Create summary data
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("7.8.2021","7.15.2021","7.22.2021","7.29.2021","8.5.2021","9.9.2021","9.15.2021"))
ggplot(data_summary) + #Fig. S9
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)

data.ecto <- data[data$fly.source=="Ecto",]
data_summary <- aggregate(percent_shared ~ as.factor(date), data.ecto,       # Create summary data
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("7.8.2021","7.15.2021","7.22.2021","7.29.2021","8.5.2021","9.9.2021","9.15.2021"))
ggplot(data_summary) + #Fig. S9
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)
			  
#### Goal: Identify ASVs that are commonly shared between fly and manure samples collected from the same facility (Fig. 6) ####
#### Goal: Plot relative abundance of commonly shared ASVs in internal/external fly vs. manure samples collected from the same facility (Fig. 7) ####

ps.Merged <- merge_samples(ps.arlington,"sample.type")
ps.Merged <- transform_sample_counts(ps.Merged, function(x) x / sum(x) )
row_names <- row.names(sample_data(ps.Merged))
sample_data(ps.Merged)$sample_id <- row_names
otus <- otu_table(ps.Merged)
otus <- as.data.frame(t(otus))
row_names <- row.names(otus)
otus$ASV <- row_names

otus <- otus[otus$ASV %in% c("54da6e3cec6a60921af7000cf5ff5618",
"aaa9219facfcb711b04e0d3c190b571a",
"4090940358ca56a1a6eaf8a8358fc385",
"e1efe25b737119a0211e4f15b35f3023",
"15d888a8d8e6c158974710b615dcbcdf",
"aa9dd0558e3ae8c191cb17006a015ac6",
"91f2efe0637f2f2062c1252c32a8e4eb",
"c27f94cb991502e7b017667c6c96f10d",
"3289d95c4a81c516faa9d263f55202d2",
"ee8bee9484867271efc76ea8c4282398",
"51dd067f061e4a29209e9a130e91fff1",
"5ce7c8bf6dc2353464c2be491f5cdf7a"), ]

otus$ASV <- factor(otus$ASV,levels = c("aa9dd0558e3ae8c191cb17006a015ac6",
"3289d95c4a81c516faa9d263f55202d2",
"ee8bee9484867271efc76ea8c4282398",
"c27f94cb991502e7b017667c6c96f10d",
"91f2efe0637f2f2062c1252c32a8e4eb",
"51dd067f061e4a29209e9a130e91fff1",
"5ce7c8bf6dc2353464c2be491f5cdf7a",
"e1efe25b737119a0211e4f15b35f3023",
"15d888a8d8e6c158974710b615dcbcdf",
"54da6e3cec6a60921af7000cf5ff5618",
"aaa9219facfcb711b04e0d3c190b571a",
"4090940358ca56a1a6eaf8a8358fc385"))

ASVNames <- c("aa9dd0558e3ae8c191cb17006a015ac6",
"3289d95c4a81c516faa9d263f55202d2",
"ee8bee9484867271efc76ea8c4282398",
"c27f94cb991502e7b017667c6c96f10d",
"91f2efe0637f2f2062c1252c32a8e4eb",
"51dd067f061e4a29209e9a130e91fff1",
"5ce7c8bf6dc2353464c2be491f5cdf7a",
"e1efe25b737119a0211e4f15b35f3023",
"15d888a8d8e6c158974710b615dcbcdf",
"54da6e3cec6a60921af7000cf5ff5618",
"aaa9219facfcb711b04e0d3c190b571a",
"4090940358ca56a1a6eaf8a8358fc385")

#nb.cols <- length(ASVNames)
#myColors <- colorRampPalette(brewer.pal(n = 8, name = "Set2"))(nb.cols)
myColors <- c("#6a51a3",
"#edf3f7",
"#b8d1e4",
"#6baed6",
"#4292c6",
"#08519c",
"#0d315d",
"#fdd0a2",
"#d94801",
"#e8e4d1",
"#d4d4d4",
"#a1a1a1")
names(myColors) <- ASVNames

ggplot(otus, aes(fill=ASV, y=Manure, x="")) + #Fig. 7
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)

ggplot(otus, aes(fill=ASV, y=Endo, x="")) + #Fig. 7
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)

ggplot(otus, aes(fill=ASV, y=Ecto, x="")) + #Fig. 7
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)

#### Goal: Plot relative abundance of commonly shared ASVs in internal/external fly vs. manure samples collected from the same facility, across different sampling dates (Fig. S11) ####

sample_data(ps)$merge_factor <- paste(sample_data(ps)$sample.type,sample_data(ps)$sampling.date, sep = "-")
ps.arlington <- subset_samples(ps, location=="DCC")
ps.Merged <- merge_samples(ps.arlington,"merge_factor")
row_names <- row.names(sample_data(ps.Merged))
sample_data(ps.Merged)$sample_id <- row_names
lst <- strsplit(row_names,'-')
v1 <- lapply(lst, `[`, 1)
v2 <- lapply(lst, `[`, 2)
sample_data(ps.Merged)$sample.type <- v1
sample_data(ps.Merged)$sampling.date <- v2

#-------------------------------------------------------#

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

#-------------------------------------------------------#

CountASVs <- function(sampleName, abundanceTable) {
  sampleAbundance <- abundanceTable %>%
    dplyr::select(sampleName) %>%
    rownames_to_column(var = "ASV")
  sampleAbundanceDF <- filter(sampleAbundance, sampleAbundance[, sampleName] > 0)
  return(sampleAbundanceDF)
}

#-------------------------------------------------------#

GetSampleSum <- function(sampleName, physeqObject) {
  named_sum <- sample_sums(physeqObject)[sampleName]
  numeric_sum <- unname(named_sum)
}

#-------------------------------------------------------#

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
      dplyr::select(as.vector("sample_id"))
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

allF1ASVTables <- allF1ASVTables[allF1ASVTables$ASV %in% c("c27f94cb991502e7b017667c6c96f10d",
"aaa9219facfcb711b04e0d3c190b571a",
"15d888a8d8e6c158974710b615dcbcdf",
"5ce7c8bf6dc2353464c2be491f5cdf7a",
"3289d95c4a81c516faa9d263f55202d2",
"51dd067f061e4a29209e9a130e91fff1",
"e1efe25b737119a0211e4f15b35f3023",
"aa9dd0558e3ae8c191cb17006a015ac6",
"4090940358ca56a1a6eaf8a8358fc385",
"91f2efe0637f2f2062c1252c32a8e4eb",
"ee8bee9484867271efc76ea8c4282398",
"54da6e3cec6a60921af7000cf5ff5618"), ]

allF1ASVTables$sample <- factor(allF1ASVTables$sample,levels = c("Manure-7.8.2021", "Endo-7.8.2021",
	"Manure-7.15.2021", "Endo-7.15.2021",
	"Manure-7.22.2021", "Endo-7.22.2021", "Ecto-7.22.2021",
	"Manure-7.29.2021", "Endo-7.29.2021", "Ecto-7.29.2021",
	"Manure-8.5.2021", "Endo-8.5.2021", "Ecto-8.5.2021",
	"Manure-8.12.2021",
	"Manure-8.19.2021",
	"Manure-8.26.2021",
	"Manure-9.2.2021",
	"Manure-9.9.2021", "Endo-9.9.2021", "Ecto-9.9.2021",
	"Manure-9.15.2021", "Endo-9.15.2021"))

ASVTax <- as.data.frame(tax_table(ps.Merged))
ASVIDs <- row.names(ASVTax)
ASVTax$ASV <- ASVIDs
ASVTax <- ASVTax[ASVTax$ASV %in% c("c27f94cb991502e7b017667c6c96f10d",
"aaa9219facfcb711b04e0d3c190b571a",
"15d888a8d8e6c158974710b615dcbcdf",
"5ce7c8bf6dc2353464c2be491f5cdf7a",
"3289d95c4a81c516faa9d263f55202d2",
"51dd067f061e4a29209e9a130e91fff1",
"e1efe25b737119a0211e4f15b35f3023",
"aa9dd0558e3ae8c191cb17006a015ac6",
"4090940358ca56a1a6eaf8a8358fc385",
"91f2efe0637f2f2062c1252c32a8e4eb",
"ee8bee9484867271efc76ea8c4282398",
"54da6e3cec6a60921af7000cf5ff5618"), ]

#	Domain	Phylum	Class	Order	Family	Genus

#aa9dd0558e3ae8c191cb17006a015ac6	Bacteria	Proteobacteria	Alphaproteobacteria	Rhizobiales	Devosiaceae	Devosia
#3289d95c4a81c516faa9d263f55202d2	Bacteria	Proteobacteria	Gammaproteobacteria	Pseudomonadales	Moraxellaceae	Acinetobacter
#ee8bee9484867271efc76ea8c4282398	Bacteria	Proteobacteria	Gammaproteobacteria	Pseudomonadales	Moraxellaceae	Acinetobacter
#c27f94cb991502e7b017667c6c96f10d	Bacteria	Proteobacteria	Gammaproteobacteria	Enterobacterales	Enterobacteriaceae	Unidentified
#91f2efe0637f2f2062c1252c32a8e4eb	Bacteria	Proteobacteria	Gammaproteobacteria	Enterobacterales	Enterobacteriaceae	Escherichia-Shigella
#51dd067f061e4a29209e9a130e91fff1	Bacteria	Proteobacteria	Gammaproteobacteria	Pseudomonadales	Pseudomonadaceae	Pseudomonas
#5ce7c8bf6dc2353464c2be491f5cdf7a	Bacteria	Proteobacteria	Gammaproteobacteria	Xanthomonadales	Xanthomonadaceae	Luteimonas
#e1efe25b737119a0211e4f15b35f3023	Bacteria	Firmicutes	Bacilli	Lactobacillales	Enterococcaceae	Enterococcus
#15d888a8d8e6c158974710b615dcbcdf	Bacteria	Firmicutes	Bacilli	Staphylococcales	Staphylococcaceae	Staphylococcus
#54da6e3cec6a60921af7000cf5ff5618	Bacteria	Actinobacteriota	Actinobacteria	Corynebacteriales	Corynebacteriaceae	Corynebacterium
#aaa9219facfcb711b04e0d3c190b571a	Bacteria	Actinobacteriota	Actinobacteria	Corynebacteriales	Dietziaceae	Dietzia
#4090940358ca56a1a6eaf8a8358fc385	Bacteria	Actinobacteriota	Actinobacteria	Micrococcales	Intrasporangiaceae	Ornithinimicrobium

#ASV5929
#ASV6721
#ASV6749
#ASV6829
#ASV6844
#ASV6409
#ASV6901
#ASV7640
#ASV7491
#ASV4727
#ASV4754
#ASV5083

allF1ASVTables$ASV <- factor(allF1ASVTables$ASV,levels = c("aa9dd0558e3ae8c191cb17006a015ac6",
"3289d95c4a81c516faa9d263f55202d2",
"ee8bee9484867271efc76ea8c4282398",
"c27f94cb991502e7b017667c6c96f10d",
"91f2efe0637f2f2062c1252c32a8e4eb",
"51dd067f061e4a29209e9a130e91fff1",
"5ce7c8bf6dc2353464c2be491f5cdf7a",
"e1efe25b737119a0211e4f15b35f3023",
"15d888a8d8e6c158974710b615dcbcdf",
"54da6e3cec6a60921af7000cf5ff5618",
"aaa9219facfcb711b04e0d3c190b571a",
"4090940358ca56a1a6eaf8a8358fc385"))

ASVNames <- c("aa9dd0558e3ae8c191cb17006a015ac6",
"3289d95c4a81c516faa9d263f55202d2",
"ee8bee9484867271efc76ea8c4282398",
"c27f94cb991502e7b017667c6c96f10d",
"91f2efe0637f2f2062c1252c32a8e4eb",
"51dd067f061e4a29209e9a130e91fff1",
"5ce7c8bf6dc2353464c2be491f5cdf7a",
"e1efe25b737119a0211e4f15b35f3023",
"15d888a8d8e6c158974710b615dcbcdf",
"54da6e3cec6a60921af7000cf5ff5618",
"aaa9219facfcb711b04e0d3c190b571a",
"4090940358ca56a1a6eaf8a8358fc385")

#nb.cols <- length(ASVNames)
#myColors <- colorRampPalette(brewer.pal(n = 8, name = "Set2"))(nb.cols)
myColors <- c("#6a51a3",
"#edf3f7",
"#b8d1e4",
"#6baed6",
"#4292c6",
"#08519c",
"#0d315d",
"#fdd0a2",
"#d94801",
"#e8e4d1",
"#d4d4d4",
"#a1a1a1")
names(myColors) <- ASVNames

F1PerspectivePlot <- ggplot(allF1ASVTables,
                                aes(x = sample, y = F1_percent, 
                                    fill = ASV)) +
  geom_bar(stat = "identity") +
scale_fill_manual(name = "ASV", values = myColors)

F1PerspectivePlot #Fig. S11
			  
#### Goal: Calculate relative abundance of manure-associated taxa that are shared with flies collected from the same facility (or not) (Fig. 5C) ####

sample_data(ps)$merge_factor <- c("Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","ArlM1-7.9.21","ArlM10-7.16.21","ArlM11-7.23.21","ArlM12-7.23.21","ArlM13-7.23.21","ArlM14-7.23.21","ArlM15-7.23.21","ArlM16-7.30.21","ArlM17-7.30.21","ArlM18-7.30.21","ArlM19-7.30.21","ArlM2-7.9.21","ArlM20-7.30.21","ArlM21-7.30.21","ArlM22-8.6.21","ArlM23-8.6.21","ArlM24-8.6.21","ArlM25-8.6.21","ArlM26-8.6.21","ArlM27-8.6.21","ArlM28-8.13.21","ArlM29-8.13.21","ArlM3-7.9.21","ArlM30-8.13.21","ArlM31-8.13.21","ArlM32-8.13.21","ArlM33-8.13.21","ArlM34-8.20.21","ArlM35-8.20.21","ArlM36-8.20.21","ArlM37-8.20.21","ArlM38-8.20.21","ArlM39-8.20.21","ArlM4-7.9.21","ArlM40-8.27.21","ArlM41-8.27.21","ArlM42-8.27.21","ArlM43-8.27.21","ArlM44-8.27.21","ArlM45-8.27.21","ArlM46-9.10.21","ArlM47-9.10.21","ArlM48-9.10.21","ArlM49-9.10.21","ArlM5-7.9.21","ArlM50-9.10.21","ArlM51-9.10.21","ArlM6-7.16.21","ArlM7-7.16.21","ArlM8-7.16.21","ArlM9-7.16.21","DCCM100-9.15.21","DCCM101-9.15.21","DCCM102-7.8.21","DCCM103-7.8.21","DCCM104-7.8.21","DCCM105-7.8.21","DCCM106-7.8.21","DCCM52-7.15.21","DCCM53-7.15.21","DCCM54-7.15.21","DCCM55-7.15.21","DCCM56-7.15.21","DCCM57-7.22.21","DCCM58-7.22.21","DCCM59-7.22.21","DCCM60-7.22.21","DCCM61-7.22.21","DCCM62-7.29.21","DCCM63-7.29.21","DCCM64-7.29.21","DCCM65-7.29.21","DCCM66-7.29.21","DCCM67-8.5.21","DCCM68-8.5.21","DCCM69-8.5.21","DCCM70-8.5.21","DCCM71-8.5.21","DCCM72-8.12.21","DCCM73-8.12.21","DCCM74-8.12.21","DCCM75-8.12.21","DCCM76-8.12.21","DCCM77-8.19.21","DCCM78-8.19.21","DCCM79-8.19.21","DCCM80-8.19.21","DCCM81-8.19.21","DCCM82-8.26.21","DCCM83-8.26.21","DCCM84-8.26.21","DCCM85-8.26.21","DCCM86-8.26.21","DCCM87-9.2.21","DCCM88-9.2.21","DCCM89-9.2.21","DCCM90-9.2.21","DCCM91-9.2.21","DCCM92-9.9.21","DCCM93-9.9.21","DCCM94-9.9.21","DCCM95-9.9.21","DCCM96-9.9.21","DCCM97-9.15.21","DCCM98-9.15.21","DCCM99-9.15.21","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly")
ps.arlington <- subset_samples(ps, location=="DCC")
ps.Merged <- merge_samples(ps.arlington,"merge_factor")
row_names <- row.names(sample_data(ps.Merged))
sample_data(ps.Merged)$sample_id <- row_names
lst <- strsplit(row_names,'-')
v1 <- lapply(lst, `[`, 1)
v2 <- lapply(lst, `[`, 2)
sample_data(ps.Merged)$pooltype <- v1
sample_data(ps.Merged)$sampling.date <- v2

#-------------------------------------------------------#

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

#-------------------------------------------------------#

CountASVs <- function(sampleName, abundanceTable) {
  sampleAbundance <- abundanceTable %>%
    dplyr::select(sampleName) %>%
    rownames_to_column(var = "ASV")
  sampleAbundanceDF <- filter(sampleAbundance, sampleAbundance[, sampleName] > 0)
  return(sampleAbundanceDF)
}

#-------------------------------------------------------#

GetSampleSum <- function(sampleName, physeqObject) {
  named_sum <- sample_sums(physeqObject)[sampleName]
  numeric_sum <- unname(named_sum)
}

#-------------------------------------------------------#

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

for (i in 1:55) {
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
      dplyr::select(as.vector("sample_id"))
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
  currentManure <- NULL

  for (i in 1:length(familyMembers)) {
      currentManure <- c(currentManure,familyMembers[[i]])
    	}

  # Put infants in a list to loop over:
  currentF1List <- list(currentManure)
  
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
	
	Fly_ASV_counts = CountASVs("Fly-Fly", abdTable)

    # merge current infant's ASV counts with mother's
    F1Table <- merge(current_F1@ASV_counts,
                         Fly_ASV_counts,
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
    
    F1Table$shared <- ifelse(!(is.na(F1Table$"Fly-Fly")), yes = TRUE, no = FALSE)
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
      dplyr::select(-c("Fly-Fly"))
    
    # Add infant sample name to table
    #infantTablePercentages <- cbind(infantTablePercentages,
                                    #"infant_sample" = currentInfantName)
    
    allF1ASVTables <- rbind(allF1ASVTables, F1TablePercentages)
  }
}

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
            file = "DCCFlyPerspectiveRAManure.csv",
            quote = FALSE, row.names = FALSE)

data=read.csv("DCCFlyPerspectiveRAManureFormatted.csv")
kruskal.test(percent_shared~as.factor(date),data=data) #NS (p = 0.4174)

mean <- mean(data$percent_shared)
ci <- 1.96*(sd(data$percent_shared) / sqrt(length(data$percent_shared)))                              
ggplot(data=NULL) + #Fig. 5C
    geom_bar( aes(x="", y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x="", ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)
				     
#### Goal: Calculate relative abundance of manure-associated taxa that are shared with flies collected from the same facility (or not), across different sampling dates (Fig. S10) ####

means <- aggregate(percent_shared ~ date, data = data, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(percent_shared ~ date, data = data, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="date")
summary$date <- factor(summary$date, levels=c("7.8.21","7.15.21","7.22.21","7.29.21","8.5.21","8.12.21","8.19.21","8.26.21","9.2.21","9.9.21","9.15.21"))
ggplot(summary, aes(x = date, y = percent_shared.x)) + #Fig. S10
    geom_errorbar(aes(ymin=percent_shared.x-percent_shared.y, ymax=percent_shared.x+percent_shared.y), width=.1) +
    geom_point()













		 
