library(phyloseq)
library(phytools)
library(ggplot2)
library(decontam)

####More pragmatic tree rooting upon import####

pick_new_outgroup <- function(tree.unrooted){
  require("magrittr")
  require("data.table")
  require("ape") # ape::Ntip
  # tablify parts of tree that we need.
  treeDT <- 
    cbind(
      data.table(tree.unrooted$edge),
      data.table(length = tree.unrooted$edge.length)
    )[1:Ntip(tree.unrooted)] %>% 
    cbind(data.table(id = tree.unrooted$tip.label))
  # Take the longest terminal branch as outgroup
  new.outgroup <- treeDT[which.max(length)]$id
  return(new.outgroup)
}

# Read in unrooted tree
unrooted_qiime_default<-read.newick(file="unrootedtree.nwk")
# the line below applies the user-defined function "pick_new_outgroup" to pick a more reasonable tree outgroup
new.outgroup = pick_new_outgroup(unrooted_qiime_default)
# the line below roots the tree using a package called "ape"
rootedTree = ape::root(unrooted_qiime_default, outgroup=new.outgroup, resolve.root=TRUE)
# the function below is the phyloseq function to make phylogenic tree
phy_tree<-phyloseq::read_tree(rootedTree, errorIfNULL=FALSE)

####Reading in files####

##OTU File
# read in the .csv file
otu_csv<-read.csv("otu_table.csv")
# set the first column with otu names to be the row names
otu_rows_fixed <- data.frame(otu_csv, row.names = 1)
# use the phyloseq function to read in the otu table making sure to tell phyloseq that each row (as opposed to each column) corresponds to a taxon
OTU_16S<-phyloseq::otu_table(otu_rows_fixed,taxa_are_rows = TRUE)

##Taxonomy File
# read in csv file
tax_csv<-read.csv("taxonomy.csv")
# set the first column with otu names to be the row names
tax_rows_fixed<-data.frame(tax_csv, row.names = 1)
# the phyloseq function for reading in a taxonomy table expects a matrix:
tax_matrix<-as.matrix(tax_rows_fixed)
# use the phyloseq function to read in a tax_table
TAX_16S<-phyloseq::tax_table(tax_matrix, errorIfNULL=TRUE)

##Metadata File
#Read in csv
meta_csv<-read.csv("Field_Samples_Metadata.csv")
# set the first column with sample names to be row names
meta_rows<-data.frame(meta_csv, row.names = 1)
# use the phyloseq function to read in sample_data
META_16S<-phyloseq::sample_data(meta_rows)

#Link together the objects to make one phyloseq object, then save it
seq_object<-phyloseq(OTU_16S,TAX_16S,META_16S,phy_tree)
#saveRDS(seq_object,"phyloseqlinkedobjects2021Field.RDS")

####Decontamination of Data####

## Plot samples 
df <- as.data.frame(sample_data(seq_object)) 
df$LibrarySize <- sample_sums(seq_object)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot2::ggplot(data=df, aes(x=Index, y=LibrarySize, color=run.type)) + geom_point() 

##Decontamination by prevalence
# Make column in metadata dataframe that identifies samples labeled as controls ("TRUE")
sample_data(seq_object)$is.neg <- sample_data(seq_object)$run.type == "Control"
contamdf.prev05 <- decontam::isContaminant(seq_object, method="prevalence", neg="is.neg", threshold=0.5) # this is a tunable threshold, using 0.5 identifies all reads that are more prevalent in controls
#table(contamdf.prev05$contaminant)
prevalence<-contamdf.prev05[which(contamdf.prev05$contaminant=='TRUE'),]
ps.noncontam <- prune_taxa(!contamdf.prev05$contaminant, seq_object) #list of all reads which meet the threshold value of 0.5 

##Decontamination by frequency
contamdf.freq <- decontam::isContaminant(ps.noncontam, method="frequency", conc="extraction.conc") ## this is also tunable with a threshold argument, using the default here
#table(contamdf.freq$contaminant)
frequency<-contamdf.freq[which(contamdf.freq$contaminant=='TRUE'),]
contaminant_reads<-c(rownames(prevalence),rownames(frequency))

# make your final phyloseq object that is nice and contaminant free!
ps.noncontam <- prune_taxa(!contamdf.freq$contaminant, ps.noncontam)

####Final Clean Ups####
greater_than_100_ASVs <- prune_samples(sample_sums(seq_object)>=100, ps.noncontam) ## get rid of any sample that contains fewer than 100 reads
no_controls<-prune_samples(sample_data(greater_than_100_ASVs)$run.type == "Sample", greater_than_100_ASVs) ## remove samples we have already used for our decontam procedure
no_controls<-subset_taxa(no_controls, Class!="Chloroplast") 
no_controls<-subset_taxa(no_controls, Family!="Mitochondria")

#saveRDS(no_controls,"ps_FieldWork2021_AJS_Final.RDS")
