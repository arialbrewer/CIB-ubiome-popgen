#Arial Brewer
#PhD Chapter 3: CIB population genetics (mitochondrial DNA)

#load packages
library(tidyverse)
library(vegan)
library(strataG)
library(adegenet)

#view strataG vignettes
browseVignettes("strataG")

#set directory
setwd("C:/Users/arial/Desktop/Ch.3 Microbiome & popgen/CIB-ubiome-popgen/Data/")

#read in data




seq.df <- dolph.strata[ c("id", "broad", "id")]
colnames(seq.df)[3] <- "D-loop"
dl.g <- df2gtypes(seq.df, ploidy = 1, sequences = dolph.seqs)
dl.g  

#label haplotypes
dl.haps <- labelHaplotypes(mito.data)
dl.haps




