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

#### Fig 7 (left, Arlington samples)

ps <- readRDS("ps_FieldWork2021_AJS_Final.rds")
ps.arlington <- subset_samples(ps, location=="Arlington")
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

ggplot(otus, aes(fill=ASV, y=Manure, x="")) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)

ggplot(otus, aes(fill=ASV, y=Endo, x="")) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)

ggplot(otus, aes(fill=ASV, y=Ecto, x="")) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)
