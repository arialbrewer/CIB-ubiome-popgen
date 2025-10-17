#Arial Brewer
#PhD Chapter 3: CIB population genetics (nuclear DNA)

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
nuc.data <- readGenData("Beluga_GTseq300_RAD_genos.csv")
  
#remove unneeded columns
nuc.data <- nuc.data %>% 
  dplyr::select(-SampleID_full,-Raw,-On_target,-OTP,-IFI) 

nuc.data <- nuc.data %>% 
  group_by(SampleID)


#create gtypes object
nuc.data <- df2gtypes(nuc.data, ploidy=2, id.col=1)

#search for duplicate genotypes
dupGenotypes(nuc.data, num.shared=0.8)


#calculate theta
theta(nuc.data)

#calculate nucleotide diversity
nucleotideDiversity(nuc.data)

#calculate heterozygosity
heterozygosity(nuc.data, type = "observed")
heterozygosity(nuc.data, type = "expected")








