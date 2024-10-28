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

x <- c(1:12)
freq <- c(11,10,9,8,7,5,5,5,5,5,5,5)

plot(x,freq)

ps <- readRDS("ps_FieldWork2021_AJS_Final.rds")
ps.arlington <- subset_samples(ps, location=="DCC")
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
