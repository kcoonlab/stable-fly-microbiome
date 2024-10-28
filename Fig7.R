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
"b40dfe8590cdabab904a9931210d2805",
"dcfd8312bf9f39a54aab3707c3b140df",
"bcaaa3d15ff438cfa14dbec5d86b0eb8",
"4090940358ca56a1a6eaf8a8358fc385",
"b77e07971c99cc0a36c99bd7cad83c16",
"aa9dd0558e3ae8c191cb17006a015ac6",
"51dd067f061e4a29209e9a130e91fff1",
"3289d95c4a81c516faa9d263f55202d2",
"350bce411d627b38ef18aeaebd7d4f69",
"ee8bee9484867271efc76ea8c4282398",
"c27f94cb991502e7b017667c6c96f10d",
"5dbd59fa8200705770f4d068d6572dc3",
"91f2efe0637f2f2062c1252c32a8e4eb",
"6207db2806b88c7f75541b8de1c42f50",
"5ce7c8bf6dc2353464c2be491f5cdf7a",
"8da7704aad448995f55464ee650a4452",
"ff27974c65664de75e6421a62a0448a3",
"15d888a8d8e6c158974710b615dcbcdf",
"945bb343f8dc3dd656b9ca15c3c0f09d",
"e0597b0224419e26dd82616c75f8c1d2",
"66402d8adb61bb82bd8d7fb58d832c10",
"e1efe25b737119a0211e4f15b35f3023",
"d56a560e43e7ab41839a1c8a10f65d36"), ]

otus$ASV <- factor(otus$ASV,levels = c("aa9dd0558e3ae8c191cb17006a015ac6",
"3289d95c4a81c516faa9d263f55202d2",
"350bce411d627b38ef18aeaebd7d4f69",
"ee8bee9484867271efc76ea8c4282398",
"6207db2806b88c7f75541b8de1c42f50",
"c27f94cb991502e7b017667c6c96f10d",
"91f2efe0637f2f2062c1252c32a8e4eb",
"5dbd59fa8200705770f4d068d6572dc3",
"51dd067f061e4a29209e9a130e91fff1",
"5ce7c8bf6dc2353464c2be491f5cdf7a",
"66402d8adb61bb82bd8d7fb58d832c10",
"e1efe25b737119a0211e4f15b35f3023",
"d56a560e43e7ab41839a1c8a10f65d36",
"8da7704aad448995f55464ee650a4452",
"ff27974c65664de75e6421a62a0448a3",
"15d888a8d8e6c158974710b615dcbcdf",
"945bb343f8dc3dd656b9ca15c3c0f09d",
"e0597b0224419e26dd82616c75f8c1d2",
"b40dfe8590cdabab904a9931210d2805",
"54da6e3cec6a60921af7000cf5ff5618",
"bcaaa3d15ff438cfa14dbec5d86b0eb8",
"aaa9219facfcb711b04e0d3c190b571a",
"b77e07971c99cc0a36c99bd7cad83c16",
"4090940358ca56a1a6eaf8a8358fc385",
"dcfd8312bf9f39a54aab3707c3b140df"))

ASVNames <- c("aa9dd0558e3ae8c191cb17006a015ac6",
"3289d95c4a81c516faa9d263f55202d2",
"350bce411d627b38ef18aeaebd7d4f69",
"ee8bee9484867271efc76ea8c4282398",
"6207db2806b88c7f75541b8de1c42f50",
"c27f94cb991502e7b017667c6c96f10d",
"91f2efe0637f2f2062c1252c32a8e4eb",
"5dbd59fa8200705770f4d068d6572dc3",
"51dd067f061e4a29209e9a130e91fff1",
"5ce7c8bf6dc2353464c2be491f5cdf7a",
"66402d8adb61bb82bd8d7fb58d832c10",
"e1efe25b737119a0211e4f15b35f3023",
"d56a560e43e7ab41839a1c8a10f65d36",
"8da7704aad448995f55464ee650a4452",
"ff27974c65664de75e6421a62a0448a3",
"15d888a8d8e6c158974710b615dcbcdf",
"945bb343f8dc3dd656b9ca15c3c0f09d",
"e0597b0224419e26dd82616c75f8c1d2",
"b40dfe8590cdabab904a9931210d2805",
"54da6e3cec6a60921af7000cf5ff5618",
"bcaaa3d15ff438cfa14dbec5d86b0eb8",
"aaa9219facfcb711b04e0d3c190b571a",
"b77e07971c99cc0a36c99bd7cad83c16",
"4090940358ca56a1a6eaf8a8358fc385",
"dcfd8312bf9f39a54aab3707c3b140df")

nb.cols <- length(ASVNames)
myColors <- colorRampPalette(brewer.pal(n = 8, name = "Set2"))(nb.cols)
names(myColors) <- ASVNames

ggplot(otus, aes(fill=ASV, y=Endo, x="")) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)

ggplot(otus, aes(fill=ASV, y=Ecto, x="")) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)

ggplot(otus, aes(fill=ASV, y=Manure, x="")) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)
