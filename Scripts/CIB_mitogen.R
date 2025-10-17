#Arial Brewer
#PhD Chapter 3: CIB population genetics (mitochondrial DNA)

#load packages
library(tidyverse)
library(vegan)
library(strataG)
library(adegenet)
library(ape)

#view strataG vignettes
browseVignettes("strataG")

#set directory
setwd("C:/Users/arial/Desktop/Ch.3 Microbiome & popgen/CIB-ubiome-popgen/Data/")

#read in data

#sample metadata
meta <- read.csv("CIB_sample_metadata.csv")

#O'Corry-Crowe haplotype data
occ.haps <- strataG::read.fasta("OCorry-Crowe_mthaps.fasta")

#short CIB mtDNA trimmed sequences
cib.short <- strataG::read.fasta("CIB_sequence_alignment_short.fasta")

#long CIB mtDNA trimmed sequences
cib.long <- strataG::read.fasta("CIB_sequence_alignment_long.fasta")

#create gtypes object
data <- df2gtypes(occ.align.haps, ploidy=1, id.col=1)



#label haplotypes
seq.df <- dolph.strata[ c("id", "broad")]
colnames(seq.df)[2] <- "D-loop"
dl.g <- df2gtypes(seq.df, ploidy = 1, sequences = occ.align.haps)
dl.g




