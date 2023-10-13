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

allF1ASVTables <- data.frame()

for (k in 1:length(validFamilies)) {
  family <- validFamilies[[k]]
  familyMembers <- family@member_list

  # Put each FamilyMember object in its own variable
  currentFly <- NULL

  for (i in 1:length(familyMembers)) {
      currentFly <- c(currentFly,familyMembers[[i]])
    	}

  # Put fly samples in a list to loop over:
  currentF1List <- list(currentFly)
  
  # Remove any NULL elements in the list (in case there is only 1 fly sample)
  currentF1List <- plyr::compact(currentF1List)

  # Loop over the fly sample list and populate the data frame:
  #   allF1ASVTables - for each fly sample, holds the ASV, 
  #     its percentage of the total reads in the sample, and whether it's shared with manure collected from the same facility
  
  for (i in 1:length(currentF1List[[1]])) {
  
    current_F1 <- currentF1List[[1]][[i]]
    #currentFlyName <- current_fly@name
    
    # build fly sample's identifier (FamGroup.1 or FamGroup.2)
    F1Designation <- current_F1@sample_designation
    F1Identifier <- paste(family@family_group, F1Designation, sep = "/")
	
	Manure_ASV_counts = CountASVs("Manure-Manure", abdTable)

    # merge current fly sample's ASV counts with manure's
    F1Table <- merge(current_F1@ASV_counts,
                         Manure_ASV_counts,
                         by = "ASV", all = TRUE)
    
    # Replace column names with the samples' sample_designations for easier
    #   plotting downstream:
    
      currentColName <- names(F1Table)[2]
      setnames(F1Table, old = currentColName, new = familyMembers[[currentColName]]@sample_designation)
    
    # Collapse table to figure out what is shared from fly sample's perspective -
    #   remove any ASVs that have a count of NA in the current sample
    F1Table <- F1Table[!(is.na(F1Table[[F1Designation]])), ]
    
    # Add a column to indicate if each ASV (row) in the collapsed table is 
    #   shared with manure or not. If the ASV has a value of NA in manure column
    #   then it cannot be shared between current fly sample and manure
    
    F1Table$shared <- ifelse(!(is.na(F1Table$"Manure-Manure")), yes = TRUE, no = FALSE)
    F1Table$sample <- F1Identifier

    # Prepare to add this table to a bigger list, but don't add until end!
    setnames(F1Table, old = F1Designation, new = "read_count")
      
    # For this fly sample, calculate percent of reads that are from ASVs that
    #   are shared with mom, and not shared with mom (keep this table).
    F1Sum <- sum(F1Table$read_count)
    F1SharedSum <- sum(F1Table$read_count[F1Table$shared == TRUE])
    F1NotSharedSum <- sum(F1Table$read_count[F1Table$shared == FALSE])
      
    # Add the flyTable (the one with ASVs as rows) to larger data frame
    F1TablePercentages <- F1Table %>%
      mutate("F1_percent" = read_count/F1Sum) %>%
      dplyr::select(-c("Manure-Manure"))
        
    allF1ASVTables <- rbind(allF1ASVTables, F1TablePercentages)
  }
}

# Group allF1ASVTables by sample and shared (T/F) then count
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
F1SharedOTUs.arlington <- F1SharedOTUs

kruskal.test(number_of_ASVs~as.factor(sample.type),data=F1SharedOTUs) #Significant (p < 0.0001)

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
F1PerspectiveRA.3$source <- ifelse(str_detect(F1PerspectiveRA.3$sample,"Endo"), "Endo", "Ecto")
F1PerspectiveRA.3 <- F1PerspectiveRA.3[F1PerspectiveRA.3$shared == "TRUE",]		 
F1SharedPercent.arlington <- F1PerspectiveRA.3
		 
kruskal.test(percent~as.factor(source),data=F1PerspectiveRA.3) #P = 0.04602
kruskal.test(percent~as.factor(date),data=F1PerspectiveRA.3) #P = 0.01786
pairwise.wilcox.test(F1PerspectiveRA.3$percent, as.factor(F1PerspectiveRA.3$date),
                 p.adjust.method = "BH", exact = FALSE) #Only significant: 7.9.2021 vs. 8.13.2021 (P = 0.030)
		 
data_summary <- aggregate(percent ~ as.factor(source), F1PerspectiveRA.3,
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
			  
F1SharedOTUs$date <- factor(F1SharedOTUs$date, levels=c("7.9.2021","7.16.2021","7.23.2021","7.30.2021","8.6.2021","8.13.2021","8.20.2021","8.27.2021","9.10.2021"))
F1SharedOTUs.Endo <- F1SharedOTUs[F1SharedOTUs$sample.type=="Endo",]
F1SharedOTUs.Ecto <- F1SharedOTUs[F1SharedOTUs$sample.type=="Ecto",]
kruskal.test(number_of_ASVs~as.factor(date),data=F1SharedOTUs.Endo) # NS (P = 0.5275)
kruskal.test(number_of_ASVs~as.factor(date),data=F1SharedOTUs.Ecto) # NS (P = 0.0781)
			  
means <- aggregate(number_of_ASVs ~ date, data = F1SharedOTUs.Ecto, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(number_of_ASVs ~ date, data = F1SharedOTUs.Ecto, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="date")
summary$date <- factor(summary$date, levels=c("7.9.2021","7.16.2021","7.23.2021","7.30.2021","8.6.2021","8.13.2021","8.20.2021","8.27.2021","9.10.2021"))
ggplot(summary, aes(x = date, y = number_of_ASVs.x)) +
    geom_errorbar(aes(ymin=number_of_ASVs.x-number_of_ASVs.y, ymax=number_of_ASVs.x+number_of_ASVs.y), width=.1) +
    geom_point()

#### Goal: Calculate the number of shared ASVs between internal/external fly samples and manure samples collected from the same facility across different sampling locations (Fig. S8) ####

means <- aggregate(number_of_ASVs ~ trap, data = F1SharedOTUs.Endo, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(number_of_ASVs ~ trap, data = F1SharedOTUs.Endo, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="trap")
summary$date <- factor(summary$trap, levels=c("B1-East", "B1-West", "B2-East", "B2-West"))
ggplot(summary, aes(x = trap, y = number_of_ASVs.x)) +
    geom_errorbar(aes(ymin=number_of_ASVs.x-number_of_ASVs.y, ymax=number_of_ASVs.x+number_of_ASVs.y), width=.1) +
    geom_point()

means <- aggregate(number_of_ASVs ~ trap, data = F1SharedOTUs.Ecto, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(number_of_ASVs ~ trap, data = F1SharedOTUs.Ecto, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="trap")
summary$date <- factor(summary$trap, levels=c("B1-East", "B1-West", "B2-East", "B2-West"))
ggplot(summary, aes(x = date, y = number_of_ASVs.x)) +
    geom_errorbar(aes(ymin=number_of_ASVs.x-number_of_ASVs.y, ymax=number_of_ASVs.x+number_of_ASVs.y), width=.1) +
    geom_point()

#### Goal: Calculate relative abundance of fly-associated taxa that are shared with manure collected from the same facility (or not), across different sampling dates (Fig. S9) ####

data.endo <- F1PerspectiveRA.3[F1PerspectiveRA.3$source=="Endo",]
data_summary <- aggregate(percent ~ as.factor(date), data.endo,
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("7.9.2021","7.16.2021","7.23.2021","7.30.2021","8.6.2021","8.13.2021","8.20.2021","8.27.2021","9.10.2021"))
ggplot(data_summary) +
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)

data.ecto <- F1PerspectiveRA.3[F1PerspectiveRA.3$source=="Ecto",]
data_summary <- aggregate(percent ~ as.factor(date), data.ecto,
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("7.9.2021","7.16.2021","7.23.2021","7.30.2021","8.6.2021","8.13.2021","8.20.2021","8.27.2021","9.10.2021"))
ggplot(data_summary) +
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)

#### Goal: Calculate relative abundance of fly-associated taxa that are shared with manure collected from the same facility (or not), across different sampling locations (Fig. S10) ####

sample_data(ps)$merge_factor <- c("Endo10A-B2-West",
"Endo10B-B2-West",
"Endo10D-B2-West",
"Endo11A-DCC",
"Endo11B-DCC",
"Endo11D-DCC",
"Endo11E-DCC",
"Endo12B-B1-East",
"Endo12C-B1-East",
"Endo12D-B1-East",
"Endo12E-B1-East",
"Endo13A-B1-West",
"Endo13B-B1-West",
"Endo13C-B1-West",
"Endo13D-B1-West",
"Endo13E-B1-West",
"Endo14A-B2-East",
"Endo14B-B2-East",
"Endo14C-B2-East",
"Endo15A-B2-West",
"Endo15B-B2-West",
"Endo15C-B2-West",
"Endo15D-B2-West",
"Endo16A-DCC",
"Endo17C-B1-East",
"Endo18A-B1-West",
"Endo18B-B1-West",
"Endo19A-B2-East",
"Endo19B-B2-East",
"Endo1A-DCC",
"Endo1B-DCC",
"Endo20A-B2-West",
"Endo20B-B2-West",
"Endo21A-DCC",
"Endo22A-B1-East",
"Endo22B-B1-East",
"Endo22C-B1-East",
"Endo22D-B1-East",
"Endo22E-B1-East",
"Endo23A-B1-West",
"Endo23B-B1-West",
"Endo23C-B1-West",
"Endo23D-B1-West",
"Endo23E-B1-West",
"Endo24A-B2-East",
"Endo24B-B2-East",
"Endo24C-B2-East",
"Endo24D-B2-East",
"Endo24E-B2-East",
"Endo25A-B2-West",
"Endo25B-B2-West",
"Endo25C-B2-West",
"Endo25D-B2-West",
"Endo27A-B1-East",
"Endo27B-B1-East",
"Endo27C-B1-East",
"Endo27D-B1-East",
"Endo27E-B1-East",
"Endo28A-B1-West",
"Endo28B-B1-West",
"Endo29A-B2-East",
"Endo29B-B2-East",
"Endo2A-B1-East",
"Endo2B-B1-East",
"Endo30A-B2-West",
"Endo30B-B2-West",
"Endo30C-B2-West",
"Endo30D-B2-West",
"Endo32A-B1-East",
"Endo32D-B1-East",
"Endo33A-B1-West",
"Endo33B-B1-West",
"Endo34A-B2-East",
"Endo34B-B2-East",
"Endo34C-B2-East",
"Endo35A-B2-West",
"Endo37A-B1-East",
"Endo37B-B1-East",
"Endo37C-B1-East",
"Endo37D-B1-East",
"Endo37E-B1-East",
"Endo38A-B1-West",
"Endo39A-B2-East",
"Endo39B-B2-East",
"Endo39C-B2-East",
"Endo39D-B2-East",
"Endo3A-B1-West",
"Endo3C-B1-West",
"Endo3D-B1-West",
"Endo3E-B1-West",
"Endo40A-B2-West",
"Endo40B-B2-West",
"Endo40C-B2-West",
"Endo41A-DCC",
"Endo42A-B1-East",
"Endo42B-B1-East",
"Endo42C-B1-East",
"Endo42D-B1-East",
"Endo42E-B1-East",
"Endo44A-B2-East",
"Endo44B-B2-East",
"Endo44C-B2-East",
"Endo44D-B2-East",
"Endo44E-B2-East",
"Endo46A-DCC",
"Endo4A-B2-East",
"Endo4B-B2-East",
"Endo4C-B2-East",
"Endo4D-B2-East",
"Endo5A-B2-West",
"Endo5B-B2-West",
"Endo6A-DCC",
"Endo6B-DCC",
"Endo6C-DCC",
"Endo6E-DCC",
"Endo7B-B1-East",
"Endo7C-B1-East",
"Endo7D-B1-East",
"Endo7E-B1-East",
"Endo8A-B1-West",
"Endo8C-B1-West",
"Endo8D-B1-West",
"Endo8E-B1-West",
"Endo9A-B2-East",
"Endo9B-B2-East",
"Endo9C-B2-East",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Manure-Manure",
"Ecto24-B2-East",
"Ecto10-B2-West",
"Ecto11-DCC",
"Ecto12-B1-East",
"Ecto13-B1-West",
"Ecto14-B2-East",
"Ecto16-DCC",
"Ecto17-B1-East",
"Ecto18-B1-West",
"Ecto19-B2-East",
"Ecto20-B2-West",
"Ecto21-DCC",
"Ecto22-B1-East",
"Ecto23-B1-West",
"Ecto25-B2-West",
"Ecto27-B1-East",
"Ecto28-B1-West",
"Ecto29-B2-East",
"Ecto30-B2-West",
"Ecto32-B1-East",
"Ecto33-B1-West",
"Ecto34-B2-East",
"Ecto35-B2-West",
"Ecto37-B1-East",
"Ecto38-B1-West",
"Ecto39-B2-East",
"Ecto40-B2-West",
"Ecto41-DCC",
"Ecto42-B1-East",
"Ecto44-B2-East",
"Ecto7-B1-East",
"Ecto9-B2-East")
ps.arlington <- subset_samples(ps, location=="Arlington")
ps.Merged <- merge_samples(ps.arlington,"merge_factor")
row_names <- row.names(sample_data(ps.Merged))
sample_data(ps.Merged)$sample_id <- row_names
lst <- strsplit(row_names,'-')
v1 <- lapply(lst, `[`, 1)
v2 <- lapply(lst, `[`, 2)
v3 <- lapply(lst, `[`, 3)
v4 <- paste(v2, v3, sep = "-")
sample_data(ps.Merged)$pooltype <- v1
sample_data(ps.Merged)$trap <- v4
sample_data(ps.Merged)$trap[sample_data(ps.Merged)$trap == "Manure-NA"] <- "Manure-Manure"

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

sampleData <- data.frame(sample_data(ps.Merged))
sampleData$sampling.date <- unlist(sampleData$trap)
sampleData$pooltype <- unlist(sampleData$pooltype)

abdTable <- as.matrix(otu_table(ps.Merged))
if (!taxa_are_rows(abdTable)) {
  abdTable <- t(abdTable)
}
abdTable <- as.data.frame(abdTable)

familyClassList <- list()

for (i in 1:140) {
  currentFamilyGroup <- sampleData[i, "trap"]
  familyNames <- names(familyClassList)
  if (!(currentFamilyGroup %in% familyNames)){
    members <- sampleData %>%
      filter(trap == currentFamilyGroup) %>%
      dplyr::select(as.vector("sample_id"))
    memberList <- as.list(members$sample_id)
    currentFamilyMembersList <- list()
    for (sample in memberList) {
      newMember <- new("FamilyMember",
                       family_group = currentFamilyGroup,
                       name = sample, 
                       sample_designation = sampleData$pooltype[sampleData$sample_id == sample],
                       sample_sum = GetSampleSum(sample, ps.Merged),
                       ASV_counts = CountASVs(sample, abdTable))
    currentFamilyMembersList[[sample]] <- newMember
    }
    totalSum <- 0
    for (member in currentFamilyMembersList) {
      sampleSum <- member@sample_sum
      totalSum <- sum(totalSum, sampleSum)
    }
    newFamily <- new("Family", 
                     family_group = currentFamilyGroup,
                     member_list = currentFamilyMembersList,
                     family_sum = totalSum)
    familyClassList[[currentFamilyGroup]] <- newFamily
  }
}
			  
allF1ASVTables <- data.frame()

for (k in 1:length(familyClassList)) {
  family <- familyClassList[[k]]
  familyMembers <- family@member_list
  currentFly <- NULL
  for (i in 1:length(familyMembers)) {
      currentFly <- c(currentFly,familyMembers[[i]])
    	}
  currentF1List <- list(currentFly)
  currentF1List <- plyr::compact(currentF1List)
  for (i in 1:length(currentF1List[[1]])) {
    current_F1 <- currentF1List[[1]][[i]]
    F1Designation <- current_F1@sample_designation
    F1Identifier <- paste(family@family_group, F1Designation, sep = "/")
	Manure_ASV_counts = CountASVs("Manure-Manure", abdTable)
    F1Table <- merge(current_F1@ASV_counts,
                         Manure_ASV_counts,
                         by = "ASV", all = TRUE)
      currentColName <- names(F1Table)[2]
      setnames(F1Table, old = currentColName, new = familyMembers[[currentColName]]@sample_designation)
    F1Table <- F1Table[!(is.na(F1Table[[F1Designation]])), ]
    F1Table$shared <- ifelse(!(is.na(F1Table$"Manure-Manure")), yes = TRUE, no = FALSE)
    F1Table$sample <- F1Identifier
    setnames(F1Table, old = F1Designation, new = "read_count")
    F1Sum <- sum(F1Table$read_count)
    F1SharedSum <- sum(F1Table$read_count[F1Table$shared == TRUE])
    F1NotSharedSum <- sum(F1Table$read_count[F1Table$shared == FALSE])
    F1TablePercentages <- F1Table %>%
      mutate("F1_percent" = read_count/F1Sum) %>%
      dplyr::select(-c("Manure-Manure"))
    allF1ASVTables <- rbind(allF1ASVTables, F1TablePercentages)
  }
}

F1PerspectiveRA <- allF1ASVTables
array <- unlist(strsplit(F1PerspectiveRA$sample, "[/]"))
trap <- array[c(TRUE, FALSE)]
F1PerspectiveRA$trap <- trap

F1PerspectiveRA.1 <- F1PerspectiveRA %>%
  group_by(trap, sample, shared) %>%
  summarise(tot_shared_reads = sum(read_count))
  
F1PerspectiveRA.2 <- F1PerspectiveRA.1 %>%
  group_by(sample) %>%
  summarise(tot_reads = sum(tot_shared_reads))

F1PerspectiveRA.3 <- merge(x = F1PerspectiveRA.1,y = F1PerspectiveRA.2,by.x = "sample",by.y = "sample", all = T)
F1PerspectiveRA.3$percent <- F1PerspectiveRA.3$tot_shared_reads/F1PerspectiveRA.3$tot_reads
F1PerspectiveRA.3$source <- ifelse(str_detect(F1PerspectiveRA.3$sample,"Endo"), "Endo", "Ecto")
F1PerspectiveRA.3 <- F1PerspectiveRA.3[F1PerspectiveRA.3$shared == "TRUE",]		 

data.endo <- F1PerspectiveRA.3[F1PerspectiveRA.3$source=="Endo",]
data_summary <- aggregate(percent ~ as.factor(trap), data.endo,
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("B1-East", "B1-West", "B2-East", "B2-West"))
ggplot(data_summary) +
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)

data.ecto <- F1PerspectiveRA.3[F1PerspectiveRA.3$source=="Ecto",]
data_summary <- aggregate(percent ~ as.factor(trap), data.ecto,
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("B1-East", "B1-West", "B2-East", "B2-West"))
ggplot(data_summary) +
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

otus$ASV <- factor(otus$ASV,levels = c("3289d95c4a81c516faa9d263f55202d2", "350bce411d627b38ef18aeaebd7d4f69","6207db2806b88c7f75541b8de1c42f50",
				      "91f2efe0637f2f2062c1252c32a8e4eb", "5dbd59fa8200705770f4d068d6572dc3", "66402d8adb61bb82bd8d7fb58d832c10",
				      "d56a560e43e7ab41839a1c8a10f65d36", "8da7704aad448995f55464ee650a4452", "ff27974c65664de75e6421a62a0448a3",
				      "15d888a8d8e6c158974710b615dcbcdf", "945bb343f8dc3dd656b9ca15c3c0f09d", "e0597b0224419e26dd82616c75f8c1d2", "b40dfe8590cdabab904a9931210d2805",
				      "bcaaa3d15ff438cfa14dbec5d86b0eb8", "aaa9219facfcb711b04e0d3c190b571a", "b77e07971c99cc0a36c99bd7cad83c16", "dcfd8312bf9f39a54aab3707c3b140df"))
				
ASVNames <- c("ASV 6644", "ASV 6664", "ASV 6775", "ASV 6767", "ASV 6760", "ASV 7558", "ASV 7573", "ASV 7392", "ASV 7403", "ASV 7414", "ASV 7420",
	     "ASV 7537", "ASV 4886", "ASV 5059", "ASV 4736", "ASV 5087", "ASV 4947")
				     
myColors <- c("#edf3f7", "#d2e2ef", "#9fc1dc", "#4292c6", "#2171b5", "#fee6ce", "#fdae6b", "#fd8d3c", "#f16913", "#d94801", "#a63603", 
              "#7f2704", "#f1efe1", "#e9e9e7", "#d4d4d4", "#bababa", "#878787")

# Manure
ggplot(otus, aes(fill=ASV, y=Manure, x="")) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)
				     
# Ecto
ggplot(otus, aes(fill=ASV, y=Ecto, x="")) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)
				     
# Endo
ggplot(otus, aes(fill=ASV, y=Ecto, x="")) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)

#### Goal: Plot relative abundance of commonly shared ASVs in internal/external fly vs. manure samples collected from the same facility, across different sampling dates (Fig. S13) ####

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

sampleData <- data.frame(sample_data(ps.Merged))
sampleData$sample.type <- unlist(sampleData$sample.type)
sampleData$sampling.date <- unlist(sampleData$sampling.date)

abdTable <- as.matrix(otu_table(ps.Merged))
if (!taxa_are_rows(abdTable)) {
  abdTable <- t(abdTable)
}
abdTable <- as.data.frame(abdTable)

familyClassList <- list()

for (i in 1:nrow(sampleData)) {
  currentFamilyGroup <- sampleData[i, "sampling.date"]
  familyNames <- names(familyClassList)
  if (!(currentFamilyGroup %in% familyNames)){
    members <- sampleData %>%
      filter(sampling.date == currentFamilyGroup) %>%
      dplyr::select(as.vector("sample_id"))
    memberList <- as.list(members$sample_id)
    currentFamilyMembersList <- list()
    for (sample in memberList) {
      newMember <- new("FamilyMember",
                       family_group = currentFamilyGroup,
                       name = sample, 
                       sample_designation = sampleData$sample.type[sampleData$sample_id == sample],
                       sample_sum = GetSampleSum(sample, ps.Merged),
                       ASV_counts = CountASVs(sample, abdTable))
    currentFamilyMembersList[[sample]] <- newMember
    }
    totalSum <- 0
    for (member in currentFamilyMembersList) {
      sampleSum <- member@sample_sum
      totalSum <- sum(totalSum, sampleSum)
    }
    newFamily <- new("Family", 
                     family_group = currentFamilyGroup,
                     member_list = currentFamilyMembersList,
                     family_sum = totalSum)
    familyClassList[[currentFamilyGroup]] <- newFamily
  }
}

#----- Determine which family groups are valid -----# 

validFamilies <- familyClassList

#----- Populate Table for Plotting Infants -----#

allF1ASVTables <- data.frame()

for (k in 1:length(validFamilies)) {
  family <- validFamilies[[k]]
  familyMembers <- family@member_list
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
  currentF1List <- list(currentManure, currentEndo, currentEcto)
  currentF1List <- plyr::compact(currentF1List)
  for (i in 1:length(currentF1List)) {
    current_F1 <- currentF1List[[i]]
    F1Identifier <- paste(current_F1@sample_designation, current_F1@family_group, sep = "-")
    F1Table <- current_F1@ASV_counts
    setnames(F1Table, old = F1Identifier, new = "read_count", skip_absent = TRUE)
    F1Sum <- sum(F1Table$read_count)
    F1TablePercentages <- F1Table %>%
      mutate("F1_percent" = read_count/F1Sum)
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
                          
allF1ASVTables$ASV <- factor(allF1ASVTables$ASV,levels = c("3289d95c4a81c516faa9d263f55202d2", "350bce411d627b38ef18aeaebd7d4f69","6207db2806b88c7f75541b8de1c42f50",
				      "91f2efe0637f2f2062c1252c32a8e4eb", "5dbd59fa8200705770f4d068d6572dc3", "66402d8adb61bb82bd8d7fb58d832c10",
				      "d56a560e43e7ab41839a1c8a10f65d36", "8da7704aad448995f55464ee650a4452", "ff27974c65664de75e6421a62a0448a3",
				      "15d888a8d8e6c158974710b615dcbcdf", "945bb343f8dc3dd656b9ca15c3c0f09d", "e0597b0224419e26dd82616c75f8c1d2", "b40dfe8590cdabab904a9931210d2805",
				      "bcaaa3d15ff438cfa14dbec5d86b0eb8", "aaa9219facfcb711b04e0d3c190b571a", "b77e07971c99cc0a36c99bd7cad83c16", "dcfd8312bf9f39a54aab3707c3b140df"))
				
ASVNames <- c("ASV 6644", "ASV 6664", "ASV 6775", "ASV 6767", "ASV 6760", "ASV 7558", "ASV 7573", "ASV 7392", "ASV 7403", "ASV 7414", "ASV 7420",
	     "ASV 7537", "ASV 4886", "ASV 5059", "ASV 4736", "ASV 5087", "ASV 4947")
				     
myColors <- c("#edf3f7", "#d2e2ef", "#9fc1dc", "#4292c6", "#2171b5", "#fee6ce", "#fdae6b", "#fd8d3c", "#f16913", "#d94801", "#a63603", 
              "#7f2704", "#f1efe1", "#e9e9e7", "#d4d4d4", "#bababa", "#878787")

F1PerspectivePlot <- ggplot(allF1ASVTables,
                                aes(x = sample, y = F1_percent, 
                                    fill = ASV)) +
  geom_bar(stat = "identity") +
scale_fill_manual(name = "ASV", values = myColors)
F1PerspectivePlot

#### Goal: Plot relative abundance of commonly shared ASVs in internal/external fly vs. manure samples collected from the same facility, across different sampling locations (Fig. S14) ####

sample_data(ps)$merge_factor <- paste(sample_data(ps)$sample.type,sample_data(ps)$trap, sep = "-")
ps.arlington <- subset_samples(ps, location=="Arlington")
ps.Merged <- merge_samples(ps.arlington,"merge_factor")
row_names <- row.names(sample_data(ps.Merged))
sample_data(ps.Merged)$sample_id <- row_names
lst <- strsplit(row_names,'-')
v1 <- lapply(lst, `[`, 1)
v2 <- lapply(lst, `[`, 2)
v3 <- lapply(lst, `[`, 3)
v4 <- paste(v2, v3, sep = "-")
sample_data(ps.Merged)$sample.type <- v1
sample_data(ps.Merged)$trap <- v4
                        
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

sampleData <- data.frame(sample_data(ps.Merged))
sampleData$sample.type <- unlist(sampleData$sample.type)
sampleData$trap <- unlist(sampleData$trap)

abdTable <- as.matrix(otu_table(ps.Merged))
if (!taxa_are_rows(abdTable)) {
  abdTable <- t(abdTable)
}
abdTable <- as.data.frame(abdTable)

familyClassList <- list()

for (i in 1:nrow(sampleData)) {
  currentFamilyGroup <- sampleData[i, "trap"]
  familyNames <- names(familyClassList)
  if (!(currentFamilyGroup %in% familyNames)){
    members <- sampleData %>%
      filter(trap == currentFamilyGroup) %>%
      dplyr::select(as.vector("sample_id"))
    memberList <- as.list(members$sample_id)
    currentFamilyMembersList <- list()
    for (sample in memberList) {
      newMember <- new("FamilyMember",
                       family_group = currentFamilyGroup,
                       name = sample, 
                       sample_designation = sampleData$sample.type[sampleData$sample_id == sample],
                       sample_sum = GetSampleSum(sample, ps.Merged),
                       ASV_counts = CountASVs(sample, abdTable))
    currentFamilyMembersList[[sample]] <- newMember
    }
    totalSum <- 0
    for (member in currentFamilyMembersList) {
      sampleSum <- member@sample_sum
      totalSum <- sum(totalSum, sampleSum)
    }
    newFamily <- new("Family", 
                     family_group = currentFamilyGroup,
                     member_list = currentFamilyMembersList,
                     family_sum = totalSum)
    familyClassList[[currentFamilyGroup]] <- newFamily
  }
}

#----- Populate Table for Plotting Infants -----#

allF1ASVTables <- data.frame()

for (k in 1:length(familyClassList)) {
  family <- familyClassList[[k]]
  familyMembers <- family@member_list
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
  currentF1List <- list(currentManure, currentEndo, currentEcto)
  currentF1List<-currentF1List[!sapply(currentF1List,is.null)]
  currentF1List <- plyr::compact(currentF1List)
  for (i in 1:length(currentF1List)) {
    current_F1 <- currentF1List[[i]]
    F1Identifier <- paste(current_F1@sample_designation, current_F1@family_group, sep = "-")
    F1Table <- current_F1@ASV_counts
    setnames(F1Table, old = F1Identifier, new = "read_count", skip_absent = TRUE)
    F1Sum <- sum(F1Table$read_count)
    F1TablePercentages <- F1Table %>%
      mutate("F1_percent" = read_count/F1Sum)
    F1TablePercentages <- cbind(F1TablePercentages,
                                    "sample" = F1Identifier)
    allF1ASVTables <- rbind(allF1ASVTables, F1TablePercentages)
  }
}

allF1ASVTables <- allF1ASVTables[allF1ASVTables$ASV %in% c("d56a560e43e7ab41839a1c8a10f65d36", "aaa9219facfcb711b04e0d3c190b571a", "3289d95c4a81c516faa9d263f55202d2", "8da7704aad448995f55464ee650a4452", "66402d8adb61bb82bd8d7fb58d832c10",
                "6207db2806b88c7f75541b8de1c42f50", "91f2efe0637f2f2062c1252c32a8e4eb", "350bce411d627b38ef18aeaebd7d4f69", "dcfd8312bf9f39a54aab3707c3b140df", "945bb343f8dc3dd656b9ca15c3c0f09d",
                "e0597b0224419e26dd82616c75f8c1d2", "bcaaa3d15ff438cfa14dbec5d86b0eb8", "5dbd59fa8200705770f4d068d6572dc3", "b77e07971c99cc0a36c99bd7cad83c16", "b40dfe8590cdabab904a9931210d2805",
                "ff27974c65664de75e6421a62a0448a3", "15d888a8d8e6c158974710b615dcbcdf"), ]

allF1ASVTables$sample <- factor(allF1ASVTables$sample,levels = c("Ecto-B1-East", "Ecto-B1-West", "Ecto-B2-East", "Ecto-B2-West", 
								 "Endo-B1-East", "Endo-B1-West", "Endo-B2-East", "Endo-B2-West",
								 "Manure-Arl-M1", "Manure-Arl-M2", "Manure-Arl-M3", "Manure-Arl-M4", "Manure-Arl-M5",
								 "Manure-Arl-sickpen"))

ASVTax <- as.data.frame(tax_table(ps.Merged))
ASVIDs <- row.names(ASVTax)
ASVTax$ASV <- ASVIDs
ASVTax <- ASVTax[ASVTax$ASV %in% c("d56a560e43e7ab41839a1c8a10f65d36", "aaa9219facfcb711b04e0d3c190b571a", "3289d95c4a81c516faa9d263f55202d2", "8da7704aad448995f55464ee650a4452", "66402d8adb61bb82bd8d7fb58d832c10",
                "6207db2806b88c7f75541b8de1c42f50", "91f2efe0637f2f2062c1252c32a8e4eb", "350bce411d627b38ef18aeaebd7d4f69", "dcfd8312bf9f39a54aab3707c3b140df", "945bb343f8dc3dd656b9ca15c3c0f09d",
                "e0597b0224419e26dd82616c75f8c1d2", "bcaaa3d15ff438cfa14dbec5d86b0eb8", "5dbd59fa8200705770f4d068d6572dc3", "b77e07971c99cc0a36c99bd7cad83c16", "b40dfe8590cdabab904a9931210d2805",
                "ff27974c65664de75e6421a62a0448a3", "15d888a8d8e6c158974710b615dcbcdf"), ]
                          
allF1ASVTables$ASV <- factor(allF1ASVTables$ASV,levels = c("3289d95c4a81c516faa9d263f55202d2", "350bce411d627b38ef18aeaebd7d4f69","6207db2806b88c7f75541b8de1c42f50",
				      "91f2efe0637f2f2062c1252c32a8e4eb", "5dbd59fa8200705770f4d068d6572dc3", "66402d8adb61bb82bd8d7fb58d832c10",
				      "d56a560e43e7ab41839a1c8a10f65d36", "8da7704aad448995f55464ee650a4452", "ff27974c65664de75e6421a62a0448a3",
				      "15d888a8d8e6c158974710b615dcbcdf", "945bb343f8dc3dd656b9ca15c3c0f09d", "e0597b0224419e26dd82616c75f8c1d2", "b40dfe8590cdabab904a9931210d2805",
				      "bcaaa3d15ff438cfa14dbec5d86b0eb8", "aaa9219facfcb711b04e0d3c190b571a", "b77e07971c99cc0a36c99bd7cad83c16", "dcfd8312bf9f39a54aab3707c3b140df"))
				
ASVNames <- c("ASV 6644", "ASV 6664", "ASV 6775", "ASV 6767", "ASV 6760", "ASV 7558", "ASV 7573", "ASV 7392", "ASV 7403", "ASV 7414", "ASV 7420",
	     "ASV 7537", "ASV 4886", "ASV 5059", "ASV 4736", "ASV 5087", "ASV 4947")
				     
myColors <- c("#edf3f7", "#d2e2ef", "#9fc1dc", "#4292c6", "#2171b5", "#fee6ce", "#fdae6b", "#fd8d3c", "#f16913", "#d94801", "#a63603", 
              "#7f2704", "#f1efe1", "#e9e9e7", "#d4d4d4", "#bababa", "#878787")

F1PerspectivePlot <- ggplot(allF1ASVTables,
                                aes(x = sample, y = F1_percent, 
                                    fill = ASV)) +
  geom_bar(stat = "identity") +
scale_fill_manual(name = "ASV", values = myColors)
F1PerspectivePlot
				     
#### Goal: Calculate relative abundance of manure-associated taxa that are shared with flies collected from the same facility (or not) (Fig. 5C) ####

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

sampleData <- data.frame(sample_data(ps.Merged))
sampleData$sampling.date <- unlist(sampleData$sampling.date)
sampleData$pooltype <- unlist(sampleData$pooltype)

abdTable <- as.matrix(otu_table(ps.Merged))
if (!taxa_are_rows(abdTable)) {
  abdTable <- t(abdTable)
}
abdTable <- as.data.frame(abdTable)

familyClassList <- list()

for (i in 1:51) {
  currentFamilyGroup <- sampleData[i, "sampling.date"]
  familyNames <- names(familyClassList)
  if (!(currentFamilyGroup %in% familyNames)){
    members <- sampleData %>%
      filter(sampling.date == currentFamilyGroup) %>%
      dplyr::select(as.vector("sample_id"))
    memberList <- as.list(members$sample_id)
    currentFamilyMembersList <- list()
    for (sample in memberList) {
      newMember <- new("FamilyMember",
                       family_group = currentFamilyGroup,
                       name = sample, 
                       sample_designation = sampleData$pooltype[sampleData$sample_id == sample],
                       sample_sum = GetSampleSum(sample, ps.Merged),
                       ASV_counts = CountASVs(sample, abdTable))
    currentFamilyMembersList[[sample]] <- newMember
    }
    totalSum <- 0
    for (member in currentFamilyMembersList) {
      sampleSum <- member@sample_sum
      totalSum <- sum(totalSum, sampleSum)
    }
    newFamily <- new("Family", 
                     family_group = currentFamilyGroup,
                     member_list = currentFamilyMembersList,
                     family_sum = totalSum)
    familyClassList[[currentFamilyGroup]] <- newFamily
  }
}

#----- Determine which family groups are valid -----# 

validFamilies <- familyClassList

#----- Populate Table for Plotting Infants -----#

allF1ASVTables <- data.frame()

for (k in 1:length(validFamilies)) {
  family <- validFamilies[[k]]
  familyMembers <- family@member_list
  currentManure <- NULL
  for (i in 1:length(familyMembers)) {
      currentManure <- c(currentManure,familyMembers[[i]])
    	}
  currentF1List <- list(currentManure)
  currentF1List <- plyr::compact(currentF1List) 
  for (i in 1:length(currentF1List[[1]])) {
    current_F1 <- currentF1List[[1]][[i]]
    F1Designation <- current_F1@sample_designation
    F1Identifier <- paste(family@family_group, F1Designation, sep = "/")
	Fly_ASV_counts = CountASVs("Fly-Fly", abdTable)
    F1Table <- merge(current_F1@ASV_counts,
                         Fly_ASV_counts,
                         by = "ASV", all = TRUE)
      currentColName <- names(F1Table)[2]
      setnames(F1Table, old = currentColName, new = familyMembers[[currentColName]]@sample_designation)
    F1Table <- F1Table[!(is.na(F1Table[[F1Designation]])), ]
    F1Table$shared <- ifelse(!(is.na(F1Table$"Fly-Fly")), yes = TRUE, no = FALSE)
    F1Table$sample <- F1Identifier
    setnames(F1Table, old = F1Designation, new = "read_count")
    F1Sum <- sum(F1Table$read_count)
    F1SharedSum <- sum(F1Table$read_count[F1Table$shared == TRUE])
    F1NotSharedSum <- sum(F1Table$read_count[F1Table$shared == FALSE])
    F1TablePercentages <- F1Table %>%
      mutate("F1_percent" = read_count/F1Sum) %>%
      dplyr::select(-c("Fly-Fly"))
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
ManureSharedPercent.arlington <- F1PerspectiveRA.3
				     
kruskal.test(percent~as.factor(date),data=F1PerspectiveRA.3) #NS (p = 1)

mean <- mean(F1PerspectiveRA.3$percent)
ci <- 1.96*(sd(F1PerspectiveRA.3$percent) / sqrt(length(F1PerspectiveRA.3$percent)))        

ggplot(data=NULL) +
    geom_bar( aes(x="", y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x="", ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)

#### Goal: Calculate relative abundance of manure-associated taxa that are shared with flies collected from the same facility (or not), across different sampling dates (Fig. S11) ####
				     
means <- aggregate(percent ~ date, data = F1PerspectiveRA.3, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(percent ~ date, data = F1PerspectiveRA.3, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="date")
summary$date <- factor(summary$date, levels=c("7.9.21","7.16.21","7.23.21","7.30.21","8.6.21","8.13.21","8.20.21","8.27.21","9.10.21"))

ggplot(summary, aes(x = date, y = percent.x)) +
    geom_errorbar(aes(ymin=percent.x-percent.y, ymax=percent.x+percent.y), width=.1) +
    geom_point()

#### Goal: Calculate relative abundance of manure-associated taxa that are shared with flies collected from the same facility (or not), across different sampling locations (Fig. S12) ####

sample_data(ps)$merge_factor <- c("Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly",
				  "ArlM1-Arl-M1","ArlM10-Arl-M5","ArlM11-Arl-M1","ArlM12-Arl-M2","ArlM13-Arl-M3","ArlM14-Arl-M4","ArlM15-Arl-M5","ArlM16-Arl-M1","ArlM17-Arl-M2","ArlM18-Arl-M3","ArlM19-Arl-M4","ArlM2-Arl-M2","ArlM20-Arl-M5","ArlM21-Arl-sickpen",
				  "ArlM22-Arl-M1","ArlM23-Arl-M2","ArlM24-Arl-M3","ArlM25-Arl-M4","ArlM26-Arl-M5","ArlM27-Arl-sickpen",
				  "ArlM28-Arl-M1","ArlM29-Arl-M2","ArlM3-Arl-M3","ArlM30-Arl-M3","ArlM31-Arl-M4","ArlM32-Arl-M5","ArlM33-Arl-sickpen",
				  "ArlM34-Arl-M1","ArlM35-Arl-M2","ArlM36-Arl-M3","ArlM37-Arl-M4","ArlM38-Arl-M5","ArlM39-Arl-sickpen",
				  "ArlM4-Arl-M4","ArlM40-Arl-M1","ArlM41-Arl-M2","ArlM42-Arl-M3","ArlM43-Arl-M4","ArlM44-Arl-M5","ArlM45-Arl-sickpen",
				  "ArlM46-Arl-M1","ArlM47-Arl-M2","ArlM48-Arl-M3","ArlM49-Arl-M4","ArlM5-Arl-M5","ArlM50-Arl-M5","ArlM51-Arl-sickpen",
				  "ArlM6-Arl-M1","ArlM7-Arl-M2","ArlM8-Arl-M3","ArlM9-Arl-M4",
				  "DCCM100-DCC-Q4","DCCM101-DCC-Outdoor",
				  "DCCM102-DCC-Q1","DCCM103-DCC-Q2","DCCM104-DCC-Q3","DCCM105-DCC-Q4","DCCM106-DCC-Outdoor",
				  "DCCM52-DCC-Q1","DCCM53-DCC-Q2","DCCM54-DCC-Q3","DCCM55-DCC-Q4","DCCM56-DCC-Outdoor",
				  "DCCM57-DCC-Q1","DCCM58-DCC-Q2","DCCM59-DCC-Q3","DCCM60-DCC-Q4","DCCM61-DCC-Outdoor",
				  "DCCM62-DCC-Q1","DCCM63-DCC-Q2","DCCM64-DCC-Q3","DCCM65-DCC-Q4","DCCM66-DCC-Outdoor",
				  "DCCM67-DCC-Q1","DCCM68-DCC-Q2","DCCM69-DCC-Q3","DCCM70-DCC-Q4","DCCM71-DCC-Outdoor",
				  "DCCM72-DCC-Q1","DCCM73-DCC-Q2","DCCM74-DCC-Q3","DCCM75-DCC-Q4","DCCM76-DCC-Outdoor",
				  "DCCM77-DCC-Q1","DCCM78-DCC-Q2","DCCM79-DCC-Q3","DCCM80-DCC-Q4","DCCM81-DCC-Outdoor",
				  "DCCM82-DCC-Q1","DCCM83-DCC-Q2","DCCM84-DCC-Q3","DCCM85-DCC-Q4","DCCM86-DCC-Outdoor",
				  "DCCM87-DCC-Q1","DCCM88-DCC-Q2","DCCM89-DCC-Q3","DCCM90-DCC-Q4","DCCM91-DCC-Outdoor",
				  "DCCM92-DCC-Q1","DCCM93-DCC-Q2","DCCM94-DCC-Q3","DCCM95-DCC-Q4","DCCM96-DCC-Outdoor",
				  "DCCM97-DCC-Q1","DCCM98-DCC-Q2","DCCM99-DCC-Q3",
				  "Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly")
ps.arlington <- subset_samples(ps, location=="Arlington")
ps.Merged <- merge_samples(ps.arlington,"merge_factor")
row_names <- row.names(sample_data(ps.Merged))
sample_data(ps.Merged)$sample_id <- row_names
lst <- strsplit(row_names,'-')
v1 <- lapply(lst, `[`, 1)
v2 <- lapply(lst, `[`, 2)
v3 <- lapply(lst, `[`, 3)
v4 <- paste(v2, v3, sep = "-")
sample_data(ps.Merged)$pooltype <- v1
sample_data(ps.Merged)$trap <- v4
sample_data(ps.Merged)$trap[sample_data(ps.Merged)$trap == "Fly-NA"] <- "Fly-Fly"

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

sampleData <- data.frame(sample_data(ps.Merged))
sampleData$trap <- unlist(sampleData$trap)
sampleData$pooltype <- unlist(sampleData$pooltype)

abdTable <- as.matrix(otu_table(ps.Merged))
if (!taxa_are_rows(abdTable)) {
  abdTable <- t(abdTable)
}
abdTable <- as.data.frame(abdTable)

familyClassList <- list()

for (i in 1:51) {
  currentFamilyGroup <- sampleData[i, "trap"]
  familyNames <- names(familyClassList)
  if (!(currentFamilyGroup %in% familyNames)){
    members <- sampleData %>%
      filter(trap == currentFamilyGroup) %>%
      dplyr::select(as.vector("sample_id"))
    memberList <- as.list(members$sample_id)
    currentFamilyMembersList <- list()
    for (sample in memberList) {
      newMember <- new("FamilyMember",
                       family_group = currentFamilyGroup,
                       name = sample, 
                       sample_designation = sampleData$pooltype[sampleData$sample_id == sample],
                       sample_sum = GetSampleSum(sample, ps.Merged),
                       ASV_counts = CountASVs(sample, abdTable))
    currentFamilyMembersList[[sample]] <- newMember
    }
    totalSum <- 0
    for (member in currentFamilyMembersList) {
      sampleSum <- member@sample_sum
      totalSum <- sum(totalSum, sampleSum)
    }
    newFamily <- new("Family", 
                     family_group = currentFamilyGroup,
                     member_list = currentFamilyMembersList,
                     family_sum = totalSum)
    familyClassList[[currentFamilyGroup]] <- newFamily
  }
}

#----- Determine which family groups are valid -----# 

validFamilies <- familyClassList

#----- Populate Table for Plotting Infants -----#

allF1ASVTables <- data.frame()

for (k in 1:length(validFamilies)) {
  family <- validFamilies[[k]]
  familyMembers <- family@member_list
  currentManure <- NULL
  for (i in 1:length(familyMembers)) {
      currentManure <- c(currentManure,familyMembers[[i]])
    	}
  currentF1List <- list(currentManure)
  currentF1List <- plyr::compact(currentF1List) 
  for (i in 1:length(currentF1List[[1]])) {
    current_F1 <- currentF1List[[1]][[i]]
    F1Designation <- current_F1@sample_designation
    F1Identifier <- paste(family@family_group, F1Designation, sep = "/")
	Fly_ASV_counts = CountASVs("Fly-Fly", abdTable)
    F1Table <- merge(current_F1@ASV_counts,
                         Fly_ASV_counts,
                         by = "ASV", all = TRUE)
      currentColName <- names(F1Table)[2]
      setnames(F1Table, old = currentColName, new = familyMembers[[currentColName]]@sample_designation)
    F1Table <- F1Table[!(is.na(F1Table[[F1Designation]])), ]
    F1Table$shared <- ifelse(!(is.na(F1Table$"Fly-Fly")), yes = TRUE, no = FALSE)
    F1Table$sample <- F1Identifier
    setnames(F1Table, old = F1Designation, new = "read_count")
    F1Sum <- sum(F1Table$read_count)
    F1SharedSum <- sum(F1Table$read_count[F1Table$shared == TRUE])
    F1NotSharedSum <- sum(F1Table$read_count[F1Table$shared == FALSE])
    F1TablePercentages <- F1Table %>%
      mutate("F1_percent" = read_count/F1Sum) %>%
      dplyr::select(-c("Fly-Fly"))
    allF1ASVTables <- rbind(allF1ASVTables, F1TablePercentages)
  }
}

F1PerspectiveRA <- allF1ASVTables
array <- unlist(strsplit(F1PerspectiveRA$sample, "[/]"))
trap <- array[c(TRUE, FALSE)]
F1PerspectiveRA$trap <- trap

F1PerspectiveRA.1 <- F1PerspectiveRA %>%
  group_by(trap, sample, shared) %>%
  summarise(tot_shared_reads = sum(read_count))
  
F1PerspectiveRA.2 <- F1PerspectiveRA.1 %>%
  group_by(sample) %>%
  summarise(tot_reads = sum(tot_shared_reads))

F1PerspectiveRA.3 <- merge(x = F1PerspectiveRA.1,y = F1PerspectiveRA.2,by.x = "sample",by.y = "sample", all = T)
F1PerspectiveRA.3$percent <- F1PerspectiveRA.3$tot_shared_reads/F1PerspectiveRA.3$tot_reads
F1PerspectiveRA.3 <- F1PerspectiveRA.3[F1PerspectiveRA.3$shared == "TRUE",]		 

data_summary <- aggregate(percent ~ as.factor(trap), F1PerspectiveRA.3,
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("Arl-M1", "Arl-M2", "Arl-M3", "Arl-M4", "Arl-M5", "Arl-sickpen"))
ggplot(data_summary) +
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)
			  
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
ps.dcc <- subset_samples(ps, location=="DCC")
ps.Merged <- merge_samples(ps.dcc,"merge_factor")
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

sampleData <- data.frame(sample_data(ps.Merged))
sampleData$sampling.date <- unlist(sampleData$sampling.date)
sampleData$pooltype <- unlist(sampleData$pooltype)

abdTable <- as.matrix(otu_table(ps.Merged))
if (!taxa_are_rows(abdTable)) {
  abdTable <- t(abdTable)
}
abdTable <- as.data.frame(abdTable)

familyClassList <- list()

for (i in 1:18) {
  currentFamilyGroup <- sampleData[i, "sampling.date"]
  familyNames <- names(familyClassList)
  if (!(currentFamilyGroup %in% familyNames)){
    members <- sampleData %>%
      filter(sampling.date == currentFamilyGroup) %>%
      dplyr::select(as.vector("sample_id"))
    memberList <- as.list(members$sample_id)
    currentFamilyMembersList <- list()
    for (sample in memberList) {
      newMember <- new("FamilyMember",
                       family_group = currentFamilyGroup,
                       name = sample, 
                       sample_designation = sampleData$pooltype[sampleData$sample_id == sample],
                       sample_sum = GetSampleSum(sample, ps.Merged),
                       ASV_counts = CountASVs(sample, abdTable))
    currentFamilyMembersList[[sample]] <- newMember
    }
    totalSum <- 0
    for (member in currentFamilyMembersList) {
      sampleSum <- member@sample_sum
      totalSum <- sum(totalSum, sampleSum)
    }
    newFamily <- new("Family", 
                     family_group = currentFamilyGroup,
                     member_list = currentFamilyMembersList,
                     family_sum = totalSum)
    familyClassList[[currentFamilyGroup]] <- newFamily
  }
}

#----- Determine which family groups are valid -----# 

validFamilies <- familyClassList

#----- Populate Table for Plotting Infants -----#

allF1ASVTables <- data.frame()
for (k in 1:length(validFamilies)) {
  family <- validFamilies[[k]]
  familyMembers <- family@member_list
  currentFly <- NULL
  for (i in 1:length(familyMembers)) {
      currentFly <- c(currentFly,familyMembers[[i]])
    	}
  currentF1List <- list(currentFly) 
  currentF1List <- plyr::compact(currentF1List)
  for (i in 1:length(currentF1List[[1]])) {
    current_F1 <- currentF1List[[1]][[i]]
    F1Designation <- current_F1@sample_designation
    F1Identifier <- paste(family@family_group, F1Designation, sep = "/")
	Manure_ASV_counts = CountASVs("Manure-Manure", abdTable)
    F1Table <- merge(current_F1@ASV_counts,
                         Manure_ASV_counts,
                         by = "ASV", all = TRUE)
      currentColName <- names(F1Table)[2]
      setnames(F1Table, old = currentColName, new = familyMembers[[currentColName]]@sample_designation)
    F1Table <- F1Table[!(is.na(F1Table[[F1Designation]])), ]
    F1Table$shared <- ifelse(!(is.na(F1Table$"Manure-Manure")), yes = TRUE, no = FALSE)
    F1Table$sample <- F1Identifier
    setnames(F1Table, old = F1Designation, new = "read_count")
    F1Sum <- sum(F1Table$read_count)
    F1SharedSum <- sum(F1Table$read_count[F1Table$shared == TRUE])
    F1NotSharedSum <- sum(F1Table$read_count[F1Table$shared == FALSE])
    F1TablePercentages <- F1Table %>%
      mutate("F1_percent" = read_count/F1Sum) %>%
      dplyr::select(-c("Manure-Manure")) 
    allF1ASVTables <- rbind(allF1ASVTables, F1TablePercentages)
  }
}

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
F1SharedOTUs.dcc <- F1SharedOTUs
			  
kruskal.test(number_of_ASVs~as.factor(sample.type),data=F1SharedOTUs) #NS (P = 0.6703)
		 
means <- aggregate(number_of_ASVs ~ sample.type, data = F1SharedOTUs, 
          FUN = function(x) c(mean = mean(x)))

cis <- aggregate(number_of_ASVs ~ sample.type, data = F1SharedOTUs, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))

summary <- merge(means,cis,by="sample.type")
summary$sample.type <- factor(summary$sample.type, levels=c("Endo", "Ecto"))

ggplot(summary, aes(x = sample.type, y = number_of_ASVs.x)) +
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
F1PerspectiveRA.3$source <- ifelse(str_detect(F1PerspectiveRA.3$sample,"Endo"), "Endo", "Ecto")
F1PerspectiveRA.3 <- F1PerspectiveRA.3[F1PerspectiveRA.3$shared == "TRUE",]
F1SharedPercent.dcc <- F1PerspectiveRA.3

kruskal.test(percent~as.factor(source),data=F1PerspectiveRA.3) #NS (P = 0.8318)
kruskal.test(percent~as.factor(date),data=F1PerspectiveRA.3) #NS (P = 0.7132)

data_summary <- aggregate(percent ~ as.factor(source), F1PerspectiveRA.3,
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

F1SharedOTUs$date <- factor(F1SharedOTUs$date, levels=c("7.8.2021","7.15.2021","7.22.2021","7.29.2021","8.5.2021","9.9.2021","9.15.2021"))
F1SharedOTUs.Endo <- F1SharedOTUs[F1SharedOTUs$sample.type=="Endo",]
F1SharedOTUs.Ecto <- F1SharedOTUs[F1SharedOTUs$sample.type=="Ecto",]
kruskal.test(number_of_ASVs~as.factor(date),data=F1SharedOTUs.Endo) #NS (P = 0.2162)
			  
means <- aggregate(number_of_ASVs ~ date, data = F1SharedOTUs.Endo, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(number_of_ASVs ~ date, data = F1SharedOTUs.Endo, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
		 
summary <- merge(means,cis,by="date")
summary$date <- factor(summary$date, levels=c("7.8.2021","7.15.2021","7.22.2021","7.29.2021","8.5.2021","9.9.2021","9.15.2021"))
		 
ggplot(summary, aes(x = date, y = number_of_ASVs.x)) +
    geom_errorbar(aes(ymin=number_of_ASVs.x-number_of_ASVs.y, ymax=number_of_ASVs.x+number_of_ASVs.y), width=.1) + 
    geom_point()
		 
ggplot(F1SharedOTUs.Ecto, aes(x=date, y=number_of_ASVs)) + geom_point()
		 
#### Goal: Calculate relative abundance of fly-associated taxa that are shared with manure collected from the same facility (or not), across different sampling dates (Fig. S9) ####

data.endo <- F1PerspectiveRA.3[F1PerspectiveRA.3$source=="Endo",]
data_summary <- aggregate(percent ~ as.factor(date), data.endo,
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("7.8.2021","7.15.2021","7.22.2021","7.29.2021","8.5.2021","9.9.2021","9.15.2021"))
ggplot(data_summary) +
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)

data.ecto <- F1PerspectiveRA.3[F1PerspectiveRA.3$source=="Ecto",]
data_summary <- aggregate(percent ~ as.factor(date), data.ecto,
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("7.8.2021","7.15.2021","7.22.2021","7.29.2021","8.5.2021","9.9.2021","9.15.2021"))
ggplot(data_summary) +
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
  filter(asvtab > 4) %>%
  arrange(asvtab)
df.asvtab.out$asv <- rownames(df.asvtab.out)
df.asvtab.out                   
tax_table(ps.Merged.manure)[df.asvtab.out$asv,]                        

ggplot(df.asvtab.out, aes(x = reorder(asv, -asvtab,sum), y = asvtab)) +
  geom_bar(fill = "#DCDCDC", stat = "identity")  
			  
#### Goal: Plot relative abundance of commonly shared ASVs in internal/external fly vs. manure samples collected from the same facility (Fig. 7) ####

ps.Merged <- merge_samples(ps.dcc,"sample.type")
ps.Merged <- transform_sample_counts(ps.Merged, function(x) x / sum(x) )
row_names <- row.names(sample_data(ps.Merged))
sample_data(ps.Merged)$sample_id <- row_names
otus <- otu_table(ps.Merged)
otus <- as.data.frame(t(otus))
row_names <- row.names(otus)
otus$ASV <- row_names

otus <- otus[otus$ASV %in% c("54da6e3cec6a60921af7000cf5ff5618","4090940358ca56a1a6eaf8a8358fc385","aa9dd0558e3ae8c191cb17006a015ac6","51dd067f061e4a29209e9a130e91fff1","ee8bee9484867271efc76ea8c4282398",
			"91f2efe0637f2f2062c1252c32a8e4eb","e1efe25b737119a0211e4f15b35f3023","3289d95c4a81c516faa9d263f55202d2","5ce7c8bf6dc2353464c2be491f5cdf7a","15d888a8d8e6c158974710b615dcbcdf",
			"aaa9219facfcb711b04e0d3c190b571a","c27f94cb991502e7b017667c6c96f10d"), ]
				     
otus$ASV <- factor(otus$ASV,levels = c("aa9dd0558e3ae8c191cb17006a015ac6", "3289d95c4a81c516faa9d263f55202d2", "ee8bee9484867271efc76ea8c4282398",
				      "c27f94cb991502e7b017667c6c96f10d", "91f2efe0637f2f2062c1252c32a8e4eb", "51dd067f061e4a29209e9a130e91fff1", "5ce7c8bf6dc2353464c2be491f5cdf7a",
				      "e1efe25b737119a0211e4f15b35f3023",
				      "15d888a8d8e6c158974710b615dcbcdf", "54da6e3cec6a60921af7000cf5ff5618",
				      "aaa9219facfcb711b04e0d3c190b571a", "4090940358ca56a1a6eaf8a8358fc385"))

ASVNames <- c("ASV 5867", "ASV 6644", "ASV 6672", "ASV 6752", "ASV 6767", "ASV 6333", "ASV 6824", "ASV 7562", "ASV 7414",
	     "ASV 4709", "ASV 4736", "ASV 5064")

myColors <- c("#6a51a3", "#edf3f7", "#b8d1e4", "#6baed6", "#4292c6", "#08519c", "#0d315d", "#fdd0a2", "#d94801", "#e8e4d1",
	"#d4d4d4", "#a1a1a1")

ggplot(otus, aes(fill=ASV, y=Manure, x="")) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)

ggplot(otus, aes(fill=ASV, y=Endo, x="")) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)

ggplot(otus, aes(fill=ASV, y=Ecto, x="")) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)

#### Goal: Plot relative abundance of commonly shared ASVs in internal/external fly vs. manure samples collected from the same facility, across different sampling dates (Fig. S15) ####

sample_data(ps)$merge_factor <- paste(sample_data(ps)$sample.type,sample_data(ps)$sampling.date, sep = "-")
ps.dcc <- subset_samples(ps, location=="DCC")
ps.Merged <- merge_samples(ps.dcc,"merge_factor")
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

sampleData <- data.frame(sample_data(ps.Merged))
sampleData$sample.type <- unlist(sampleData$sample.type)
sampleData$sampling.date <- unlist(sampleData$sampling.date)

abdTable <- as.matrix(otu_table(ps.Merged))
if (!taxa_are_rows(abdTable)) {
  abdTable <- t(abdTable)
}
abdTable <- as.data.frame(abdTable)

familyClassList <- list()

for (i in 1:nrow(sampleData)) {
  currentFamilyGroup <- sampleData[i, "sampling.date"]
  familyNames <- names(familyClassList)
  if (!(currentFamilyGroup %in% familyNames)){
    members <- sampleData %>%
      filter(sampling.date == currentFamilyGroup) %>%
      dplyr::select(as.vector("sample_id"))
    memberList <- as.list(members$sample_id)
    currentFamilyMembersList <- list()
    for (sample in memberList) {
      newMember <- new("FamilyMember",
                       family_group = currentFamilyGroup,
                       name = sample, 
                       sample_designation = sampleData$sample.type[sampleData$sample_id == sample],
                       sample_sum = GetSampleSum(sample, ps.Merged),
                       ASV_counts = CountASVs(sample, abdTable))
    currentFamilyMembersList[[sample]] <- newMember
    }
    totalSum <- 0
    for (member in currentFamilyMembersList) {
      sampleSum <- member@sample_sum
      totalSum <- sum(totalSum, sampleSum)
    }
    newFamily <- new("Family", 
                     family_group = currentFamilyGroup,
                     member_list = currentFamilyMembersList,
                     family_sum = totalSum)
    familyClassList[[currentFamilyGroup]] <- newFamily
  }
}

#----- Determine which family groups are valid -----# 

validFamilies <- familyClassList

allF1ASVTables <- data.frame()

for (k in 1:length(validFamilies)) {
  family <- validFamilies[[k]]
  familyMembers <- family@member_list
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
  currentF1List <- list(currentManure, currentEndo, currentEcto)
  currentF1List <- plyr::compact(currentF1List)
  for (i in 1:length(currentF1List)) {    
    current_F1 <- currentF1List[[i]]
    F1Identifier <- paste(current_F1@sample_designation, current_F1@family_group, sep = "-")
    F1Table <- current_F1@ASV_counts
    setnames(F1Table, old = F1Identifier, new = "read_count", skip_absent = TRUE)
    F1Sum <- sum(F1Table$read_count)
    F1TablePercentages <- F1Table %>%
      mutate("F1_percent" = read_count/F1Sum)
    F1TablePercentages <- cbind(F1TablePercentages,
                                    "sample" = F1Identifier)
    allF1ASVTables <- rbind(allF1ASVTables, F1TablePercentages)
  }
}

allF1ASVTables <- allF1ASVTables[allF1ASVTables$ASV %in% c("54da6e3cec6a60921af7000cf5ff5618","4090940358ca56a1a6eaf8a8358fc385","aa9dd0558e3ae8c191cb17006a015ac6","51dd067f061e4a29209e9a130e91fff1","ee8bee9484867271efc76ea8c4282398",
			"91f2efe0637f2f2062c1252c32a8e4eb","e1efe25b737119a0211e4f15b35f3023","3289d95c4a81c516faa9d263f55202d2","5ce7c8bf6dc2353464c2be491f5cdf7a","15d888a8d8e6c158974710b615dcbcdf",
			"aaa9219facfcb711b04e0d3c190b571a","c27f94cb991502e7b017667c6c96f10d"), ]

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
ASVTax <- ASVTax[ASVTax$ASV %in% c("54da6e3cec6a60921af7000cf5ff5618","4090940358ca56a1a6eaf8a8358fc385","aa9dd0558e3ae8c191cb17006a015ac6","51dd067f061e4a29209e9a130e91fff1","ee8bee9484867271efc76ea8c4282398",
			"91f2efe0637f2f2062c1252c32a8e4eb","e1efe25b737119a0211e4f15b35f3023","3289d95c4a81c516faa9d263f55202d2","5ce7c8bf6dc2353464c2be491f5cdf7a","15d888a8d8e6c158974710b615dcbcdf",
			"aaa9219facfcb711b04e0d3c190b571a","c27f94cb991502e7b017667c6c96f10d"), ]

allF1ASVTables$ASV <- factor(allF1ASVTables$ASV,levels = c("aa9dd0558e3ae8c191cb17006a015ac6", "3289d95c4a81c516faa9d263f55202d2", "ee8bee9484867271efc76ea8c4282398",
				      "c27f94cb991502e7b017667c6c96f10d", "91f2efe0637f2f2062c1252c32a8e4eb", "51dd067f061e4a29209e9a130e91fff1", "5ce7c8bf6dc2353464c2be491f5cdf7a",
				      "e1efe25b737119a0211e4f15b35f3023",
				      "15d888a8d8e6c158974710b615dcbcdf", "54da6e3cec6a60921af7000cf5ff5618",
				      "aaa9219facfcb711b04e0d3c190b571a", "4090940358ca56a1a6eaf8a8358fc385"))

ASVNames <- c("ASV 5867", "ASV 6644", "ASV 6672", "ASV 6752", "ASV 6767", "ASV 6333", "ASV 6824", "ASV 7562", "ASV 7414",
	     "ASV 4709", "ASV 4736", "ASV 5064")

myColors <- c("#6a51a3", "#edf3f7", "#b8d1e4", "#6baed6", "#4292c6", "#08519c", "#0d315d", "#fdd0a2", "#d94801", "#e8e4d1",
	"#d4d4d4", "#a1a1a1")

F1PerspectivePlot <- ggplot(allF1ASVTables,
                                aes(x = sample, y = F1_percent, 
                                    fill = ASV)) +
  geom_bar(stat = "identity") +
scale_fill_manual(name = "ASV", values = myColors)

F1PerspectivePlot

#### Goal: Plot relative abundance of commonly shared ASVs in internal/external fly vs. manure samples collected from the same facility, across different sampling locations (Fig. S16) ####
				     
sample_data(ps)$merge_factor <- paste(sample_data(ps)$sample.type,sample_data(ps)$trap, sep = "-")
ps.dcc <- subset_samples(ps, location=="DCC" & sample.type=="Manure")
ps.Merged <- merge_samples(ps.dcc,"merge_factor")
row_names <- row.names(sample_data(ps.Merged))
sample_data(ps.Merged)$sample_id <- row_names
lst <- strsplit(row_names,'-')
v1 <- lapply(lst, `[`, 1)
v2 <- lapply(lst, `[`, 2)
v3 <- lapply(lst, `[`, 3)
v4 <- paste(v2, v3, sep = "-")
sample_data(ps.Merged)$sample.type <- v1
sample_data(ps.Merged)$trap <- v4
                        
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

sampleData <- data.frame(sample_data(ps.Merged))
sampleData$sample.type <- unlist(sampleData$sample.type)
sampleData$trap <- unlist(sampleData$trap)

abdTable <- as.matrix(otu_table(ps.Merged))
if (!taxa_are_rows(abdTable)) {
  abdTable <- t(abdTable)
}
abdTable <- as.data.frame(abdTable)

familyClassList <- list()

for (i in 1:nrow(sampleData)) {
  currentFamilyGroup <- sampleData[i, "trap"]
  familyNames <- names(familyClassList)
  if (!(currentFamilyGroup %in% familyNames)){
    members <- sampleData %>%
      filter(trap == currentFamilyGroup) %>%
      dplyr::select(as.vector("sample_id"))
    memberList <- as.list(members$sample_id)
    currentFamilyMembersList <- list()
    for (sample in memberList) {
      newMember <- new("FamilyMember",
                       family_group = currentFamilyGroup,
                       name = sample, 
                       sample_designation = sampleData$sample.type[sampleData$sample_id == sample],
                       sample_sum = GetSampleSum(sample, ps.Merged),
                       ASV_counts = CountASVs(sample, abdTable))
    currentFamilyMembersList[[sample]] <- newMember
    }
    totalSum <- 0
    for (member in currentFamilyMembersList) {
      sampleSum <- member@sample_sum
      totalSum <- sum(totalSum, sampleSum)
    }
    newFamily <- new("Family", 
                     family_group = currentFamilyGroup,
                     member_list = currentFamilyMembersList,
                     family_sum = totalSum)
    familyClassList[[currentFamilyGroup]] <- newFamily
  }
}

#----- Populate Table for Plotting Infants -----#

allF1ASVTables <- data.frame()

for (k in 1:length(familyClassList)) {
  family <- familyClassList[[k]]
  familyMembers <- family@member_list
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
  currentF1List <- list(currentManure, currentEndo, currentEcto)
  currentF1List<-currentF1List[!sapply(currentF1List,is.null)]
  currentF1List <- plyr::compact(currentF1List)
  for (i in 1:length(currentF1List)) {
    current_F1 <- currentF1List[[i]]
    F1Identifier <- paste(current_F1@sample_designation, current_F1@family_group, sep = "-")
    F1Table <- current_F1@ASV_counts
    setnames(F1Table, old = F1Identifier, new = "read_count", skip_absent = TRUE)
    F1Sum <- sum(F1Table$read_count)
    F1TablePercentages <- F1Table %>%
      mutate("F1_percent" = read_count/F1Sum)
    F1TablePercentages <- cbind(F1TablePercentages,
                                    "sample" = F1Identifier)
    allF1ASVTables <- rbind(allF1ASVTables, F1TablePercentages)
  }
}

allF1ASVTables <- allF1ASVTables[allF1ASVTables$ASV %in% c("54da6e3cec6a60921af7000cf5ff5618","4090940358ca56a1a6eaf8a8358fc385","aa9dd0558e3ae8c191cb17006a015ac6","51dd067f061e4a29209e9a130e91fff1","ee8bee9484867271efc76ea8c4282398",
			"91f2efe0637f2f2062c1252c32a8e4eb","e1efe25b737119a0211e4f15b35f3023","3289d95c4a81c516faa9d263f55202d2","5ce7c8bf6dc2353464c2be491f5cdf7a","15d888a8d8e6c158974710b615dcbcdf",
			"aaa9219facfcb711b04e0d3c190b571a","c27f94cb991502e7b017667c6c96f10d"), ]

allF1ASVTables$sample <- factor(allF1ASVTables$sample,levels = c("Manure-DCC-Outdoor", "Manure-DCC-Q1", "Manure-DCC-Q2", "Manure-DCC-Q3", "Manure-DCC-Q4"))

ASVTax <- as.data.frame(tax_table(ps.Merged))
ASVIDs <- row.names(ASVTax)
ASVTax$ASV <- ASVIDs
ASVTax <- ASVTax[ASVTax$ASV %in% c("54da6e3cec6a60921af7000cf5ff5618","4090940358ca56a1a6eaf8a8358fc385","aa9dd0558e3ae8c191cb17006a015ac6","51dd067f061e4a29209e9a130e91fff1","ee8bee9484867271efc76ea8c4282398",
			"91f2efe0637f2f2062c1252c32a8e4eb","e1efe25b737119a0211e4f15b35f3023","3289d95c4a81c516faa9d263f55202d2","5ce7c8bf6dc2353464c2be491f5cdf7a","15d888a8d8e6c158974710b615dcbcdf",
			"aaa9219facfcb711b04e0d3c190b571a","c27f94cb991502e7b017667c6c96f10d"), ]
                          
allF1ASVTables$ASV <- factor(allF1ASVTables$ASV,levels = c("aa9dd0558e3ae8c191cb17006a015ac6", "3289d95c4a81c516faa9d263f55202d2", "ee8bee9484867271efc76ea8c4282398",
				      "c27f94cb991502e7b017667c6c96f10d", "91f2efe0637f2f2062c1252c32a8e4eb", "51dd067f061e4a29209e9a130e91fff1", "5ce7c8bf6dc2353464c2be491f5cdf7a",
				      "e1efe25b737119a0211e4f15b35f3023",
				      "15d888a8d8e6c158974710b615dcbcdf", "54da6e3cec6a60921af7000cf5ff5618",
				      "aaa9219facfcb711b04e0d3c190b571a", "4090940358ca56a1a6eaf8a8358fc385"))

ASVNames <- c("ASV 5867", "ASV 6644", "ASV 6672", "ASV 6752", "ASV 6767", "ASV 6333", "ASV 6824", "ASV 7562", "ASV 7414",
	     "ASV 4709", "ASV 4736", "ASV 5064")

myColors <- c("#6a51a3", "#edf3f7", "#b8d1e4", "#6baed6", "#4292c6", "#08519c", "#0d315d", "#fdd0a2", "#d94801", "#e8e4d1",
	"#d4d4d4", "#a1a1a1")

F1PerspectivePlot <- ggplot(allF1ASVTables,
                                aes(x = sample, y = F1_percent, 
                                    fill = ASV)) +
  geom_bar(stat = "identity") +
scale_fill_manual(name = "ASV", values = myColors)
F1PerspectivePlot
				     
#### Goal: Calculate relative abundance of manure-associated taxa that are shared with flies collected from the same facility (or not) (Fig. 5C) ####

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

sampleData <- data.frame(sample_data(ps.Merged))
sampleData$sampling.date <- unlist(sampleData$sampling.date)
sampleData$pooltype <- unlist(sampleData$pooltype)
abdTable <- as.matrix(otu_table(ps.Merged))
if (!taxa_are_rows(abdTable)) {
  abdTable <- t(abdTable)
}
abdTable <- as.data.frame(abdTable)

familyClassList <- list()

for (i in 1:55) {
  currentFamilyGroup <- sampleData[i, "sampling.date"]
  familyNames <- names(familyClassList)
  if (!(currentFamilyGroup %in% familyNames)){
    members <- sampleData %>%
      filter(sampling.date == currentFamilyGroup) %>%
      dplyr::select(as.vector("sample_id"))
    memberList <- as.list(members$sample_id)
    currentFamilyMembersList <- list()
    for (sample in memberList) {
      newMember <- new("FamilyMember",
                       family_group = currentFamilyGroup,
                       name = sample, 
                       sample_designation = sampleData$pooltype[sampleData$sample_id == sample],
                       sample_sum = GetSampleSum(sample, ps.Merged),
                       ASV_counts = CountASVs(sample, abdTable))
    currentFamilyMembersList[[sample]] <- newMember
    }
    totalSum <- 0
    for (member in currentFamilyMembersList) {
      sampleSum <- member@sample_sum
      totalSum <- sum(totalSum, sampleSum)
    }
    newFamily <- new("Family", 
                     family_group = currentFamilyGroup,
                     member_list = currentFamilyMembersList,
                     family_sum = totalSum)
    familyClassList[[currentFamilyGroup]] <- newFamily
  }
}

#----- Determine which family groups are valid -----# 

validFamilies <- familyClassList

#----- Populate Table for Plotting Infants -----#

allF1ASVTables <- data.frame()

for (k in 1:length(validFamilies)) {
  family <- validFamilies[[k]]
  familyMembers <- family@member_list
  currentManure <- NULL
  for (i in 1:length(familyMembers)) {
      currentManure <- c(currentManure,familyMembers[[i]])
    	}
  currentF1List <- list(currentManure)
  currentF1List <- plyr::compact(currentF1List)
  for (i in 1:length(currentF1List[[1]])) {
    current_F1 <- currentF1List[[1]][[i]]
    F1Designation <- current_F1@sample_designation
    F1Identifier <- paste(family@family_group, F1Designation, sep = "/")
	Fly_ASV_counts = CountASVs("Fly-Fly", abdTable)
    F1Table <- merge(current_F1@ASV_counts,
                         Fly_ASV_counts,
                         by = "ASV", all = TRUE)
      currentColName <- names(F1Table)[2]
      setnames(F1Table, old = currentColName, new = familyMembers[[currentColName]]@sample_designation)
    F1Table <- F1Table[!(is.na(F1Table[[F1Designation]])), ]
    F1Table$shared <- ifelse(!(is.na(F1Table$"Fly-Fly")), yes = TRUE, no = FALSE)
    F1Table$sample <- F1Identifier
    setnames(F1Table, old = F1Designation, new = "read_count")
    F1Sum <- sum(F1Table$read_count)
    F1SharedSum <- sum(F1Table$read_count[F1Table$shared == TRUE])
    F1NotSharedSum <- sum(F1Table$read_count[F1Table$shared == FALSE])
    F1TablePercentages <- F1Table %>%
      mutate("F1_percent" = read_count/F1Sum) %>%
      dplyr::select(-c("Fly-Fly")) 
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
ManureSharedPercent.dcc <- F1PerspectiveRA.3
				     
kruskal.test(percent~as.factor(date),data=F1PerspectiveRA.3) #NS (P = 0.4213)

mean <- mean(F1PerspectiveRA.3$percent)
ci <- 1.96*(sd(F1PerspectiveRA.3$percent) / sqrt(length(F1PerspectiveRA.3$percent)))                              
ggplot(data=NULL) +
    geom_bar( aes(x="", y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x="", ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)
				     
#### Goal: Calculate relative abundance of manure-associated taxa that are shared with flies collected from the same facility (or not), across different sampling dates (Fig. S11) ####

means <- aggregate(percent ~ date, data = F1PerspectiveRA.3, 
          FUN = function(x) c(mean = mean(x)))
cis <- aggregate(percent ~ date, data = F1PerspectiveRA.3, 
          FUN = function(x) c(ci = 1.96*(sd(x) / sqrt(length(x)))))
summary <- merge(means,cis,by="date")
summary$date <- factor(summary$date, levels=c("7.8.21","7.15.21","7.22.21","7.29.21","8.5.21","8.12.21","8.19.21","8.26.21","9.2.21","9.9.21","9.15.21"))
ggplot(summary, aes(x = date, y = percent.x)) + #Fig. S10
    geom_errorbar(aes(ymin=percent.x-percent.y, ymax=percent.x+percent.y), width=.1) +
    geom_point()

#### Goal: Calculate relative abundance of manure-associated taxa that are shared with flies collected from the same facility (or not), across different sampling locations (Fig. S12) ####

sample_data(ps)$merge_factor <- c("Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly",
				  "ArlM1-Arl-M1","ArlM10-Arl-M5","ArlM11-Arl-M1","ArlM12-Arl-M2","ArlM13-Arl-M3","ArlM14-Arl-M4","ArlM15-Arl-M5","ArlM16-Arl-M1","ArlM17-Arl-M2","ArlM18-Arl-M3","ArlM19-Arl-M4","ArlM2-Arl-M2","ArlM20-Arl-M5","ArlM21-Arl-sickpen",
				  "ArlM22-Arl-M1","ArlM23-Arl-M2","ArlM24-Arl-M3","ArlM25-Arl-M4","ArlM26-Arl-M5","ArlM27-Arl-sickpen",
				  "ArlM28-Arl-M1","ArlM29-Arl-M2","ArlM3-Arl-M3","ArlM30-Arl-M3","ArlM31-Arl-M4","ArlM32-Arl-M5","ArlM33-Arl-sickpen",
				  "ArlM34-Arl-M1","ArlM35-Arl-M2","ArlM36-Arl-M3","ArlM37-Arl-M4","ArlM38-Arl-M5","ArlM39-Arl-sickpen",
				  "ArlM4-Arl-M4","ArlM40-Arl-M1","ArlM41-Arl-M2","ArlM42-Arl-M3","ArlM43-Arl-M4","ArlM44-Arl-M5","ArlM45-Arl-sickpen",
				  "ArlM46-Arl-M1","ArlM47-Arl-M2","ArlM48-Arl-M3","ArlM49-Arl-M4","ArlM5-Arl-M5","ArlM50-Arl-M5","ArlM51-Arl-sickpen",
				  "ArlM6-Arl-M1","ArlM7-Arl-M2","ArlM8-Arl-M3","ArlM9-Arl-M4",
				  "DCCM100-DCC-Q4","DCCM101-DCC-Outdoor",
				  "DCCM102-DCC-Q1","DCCM103-DCC-Q2","DCCM104-DCC-Q3","DCCM105-DCC-Q4","DCCM106-DCC-Outdoor",
				  "DCCM52-DCC-Q1","DCCM53-DCC-Q2","DCCM54-DCC-Q3","DCCM55-DCC-Q4","DCCM56-DCC-Outdoor",
				  "DCCM57-DCC-Q1","DCCM58-DCC-Q2","DCCM59-DCC-Q3","DCCM60-DCC-Q4","DCCM61-DCC-Outdoor",
				  "DCCM62-DCC-Q1","DCCM63-DCC-Q2","DCCM64-DCC-Q3","DCCM65-DCC-Q4","DCCM66-DCC-Outdoor",
				  "DCCM67-DCC-Q1","DCCM68-DCC-Q2","DCCM69-DCC-Q3","DCCM70-DCC-Q4","DCCM71-DCC-Outdoor",
				  "DCCM72-DCC-Q1","DCCM73-DCC-Q2","DCCM74-DCC-Q3","DCCM75-DCC-Q4","DCCM76-DCC-Outdoor",
				  "DCCM77-DCC-Q1","DCCM78-DCC-Q2","DCCM79-DCC-Q3","DCCM80-DCC-Q4","DCCM81-DCC-Outdoor",
				  "DCCM82-DCC-Q1","DCCM83-DCC-Q2","DCCM84-DCC-Q3","DCCM85-DCC-Q4","DCCM86-DCC-Outdoor",
				  "DCCM87-DCC-Q1","DCCM88-DCC-Q2","DCCM89-DCC-Q3","DCCM90-DCC-Q4","DCCM91-DCC-Outdoor",
				  "DCCM92-DCC-Q1","DCCM93-DCC-Q2","DCCM94-DCC-Q3","DCCM95-DCC-Q4","DCCM96-DCC-Outdoor",
				  "DCCM97-DCC-Q1","DCCM98-DCC-Q2","DCCM99-DCC-Q3",
				  "Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly","Fly-Fly")
ps.dcc <- subset_samples(ps, location=="DCC")
ps.Merged <- merge_samples(ps.dcc,"merge_factor")
row_names <- row.names(sample_data(ps.Merged))
sample_data(ps.Merged)$sample_id <- row_names
lst <- strsplit(row_names,'-')
v1 <- lapply(lst, `[`, 1)
v2 <- lapply(lst, `[`, 2)
v3 <- lapply(lst, `[`, 3)
v4 <- paste(v2, v3, sep = "-")
sample_data(ps.Merged)$pooltype <- v1
sample_data(ps.Merged)$trap <- v4
sample_data(ps.Merged)$trap[sample_data(ps.Merged)$trap == "Fly-NA"] <- "Fly-Fly"

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

sampleData <- data.frame(sample_data(ps.Merged))
sampleData$trap <- unlist(sampleData$trap)
sampleData$pooltype <- unlist(sampleData$pooltype)

abdTable <- as.matrix(otu_table(ps.Merged))
if (!taxa_are_rows(abdTable)) {
  abdTable <- t(abdTable)
}
abdTable <- as.data.frame(abdTable)

familyClassList <- list()

for (i in 1:51) {
  currentFamilyGroup <- sampleData[i, "trap"]
  familyNames <- names(familyClassList)
  if (!(currentFamilyGroup %in% familyNames)){
    members <- sampleData %>%
      filter(trap == currentFamilyGroup) %>%
      dplyr::select(as.vector("sample_id"))
    memberList <- as.list(members$sample_id)
    currentFamilyMembersList <- list()
    for (sample in memberList) {
      newMember <- new("FamilyMember",
                       family_group = currentFamilyGroup,
                       name = sample, 
                       sample_designation = sampleData$pooltype[sampleData$sample_id == sample],
                       sample_sum = GetSampleSum(sample, ps.Merged),
                       ASV_counts = CountASVs(sample, abdTable))
    currentFamilyMembersList[[sample]] <- newMember
    }
    totalSum <- 0
    for (member in currentFamilyMembersList) {
      sampleSum <- member@sample_sum
      totalSum <- sum(totalSum, sampleSum)
    }
    newFamily <- new("Family", 
                     family_group = currentFamilyGroup,
                     member_list = currentFamilyMembersList,
                     family_sum = totalSum)
    familyClassList[[currentFamilyGroup]] <- newFamily
  }
}

#----- Determine which family groups are valid -----# 

validFamilies <- familyClassList

#----- Populate Table for Plotting Infants -----#

allF1ASVTables <- data.frame()

for (k in 1:length(validFamilies)) {
  family <- validFamilies[[k]]
  familyMembers <- family@member_list
  currentManure <- NULL
  for (i in 1:length(familyMembers)) {
      currentManure <- c(currentManure,familyMembers[[i]])
    	}
  currentF1List <- list(currentManure)
  currentF1List <- plyr::compact(currentF1List) 
  for (i in 1:length(currentF1List[[1]])) {
    current_F1 <- currentF1List[[1]][[i]]
    F1Designation <- current_F1@sample_designation
    F1Identifier <- paste(family@family_group, F1Designation, sep = "/")
	Fly_ASV_counts = CountASVs("Fly-Fly", abdTable)
    F1Table <- merge(current_F1@ASV_counts,
                         Fly_ASV_counts,
                         by = "ASV", all = TRUE)
      currentColName <- names(F1Table)[2]
      setnames(F1Table, old = currentColName, new = familyMembers[[currentColName]]@sample_designation)
    F1Table <- F1Table[!(is.na(F1Table[[F1Designation]])), ]
    F1Table$shared <- ifelse(!(is.na(F1Table$"Fly-Fly")), yes = TRUE, no = FALSE)
    F1Table$sample <- F1Identifier
    setnames(F1Table, old = F1Designation, new = "read_count")
    F1Sum <- sum(F1Table$read_count)
    F1SharedSum <- sum(F1Table$read_count[F1Table$shared == TRUE])
    F1NotSharedSum <- sum(F1Table$read_count[F1Table$shared == FALSE])
    F1TablePercentages <- F1Table %>%
      mutate("F1_percent" = read_count/F1Sum) %>%
      dplyr::select(-c("Fly-Fly"))
    allF1ASVTables <- rbind(allF1ASVTables, F1TablePercentages)
  }
}

F1PerspectiveRA <- allF1ASVTables
array <- unlist(strsplit(F1PerspectiveRA$sample, "[/]"))
trap <- array[c(TRUE, FALSE)]
F1PerspectiveRA$trap <- trap

F1PerspectiveRA.1 <- F1PerspectiveRA %>%
  group_by(trap, sample, shared) %>%
  summarise(tot_shared_reads = sum(read_count))
  
F1PerspectiveRA.2 <- F1PerspectiveRA.1 %>%
  group_by(sample) %>%
  summarise(tot_reads = sum(tot_shared_reads))

F1PerspectiveRA.3 <- merge(x = F1PerspectiveRA.1,y = F1PerspectiveRA.2,by.x = "sample",by.y = "sample", all = T)
F1PerspectiveRA.3$percent <- F1PerspectiveRA.3$tot_shared_reads/F1PerspectiveRA.3$tot_reads
F1PerspectiveRA.3 <- F1PerspectiveRA.3[F1PerspectiveRA.3$shared == "TRUE",]		 

data_summary <- aggregate(percent ~ as.factor(trap), F1PerspectiveRA.3,
                          function(x) c(mean = mean(x),
                                        se = sd(x) / sqrt(length(x))))
data_summary <- data.frame(group = data_summary[,1], data_summary[,2])
data_summary$ci <- 1.96*data_summary$se
data_summary$group <- factor(data_summary$group, levels=c("DCC-Q1", "DCC-Q2", "DCC-Q3", "DCC-Q4", "DCC-Outdoor"))
ggplot(data_summary) +
    geom_bar( aes(x=group, y=mean), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean-ci, ymax=mean+ci), width=0.4, colour="orange", alpha=0.9, size=1.3) +
    ylim(0,1)
