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

otus <- otus[otus$ASV %in% c("15d888a8d8e6c158974710b615dcbcdf",
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

otus$ASV <- factor(otus$ASV,levels = c("3289d95c4a81c516faa9d263f55202d2",
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

ggplot(otus, aes(fill=ASV, y=Endo, x="")) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)

ggplot(otus, aes(fill=ASV, y=Ecto, x="")) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)

ggplot(otus, aes(fill=ASV, y=Manure, x="")) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)                                    

#### Fig 7 (right, DCC samples)

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

ggplot(otus, aes(fill=ASV, y=Manure, x="")) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(name = "ASV", values = myColors)
                                     
