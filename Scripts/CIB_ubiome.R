#Arial Brewer
#PhD Chapter 3: CIB epidermal microbiome

#load packages
library(dada2)
library(phyloseq)
library(tidyverse)

######################DADA2 pipeline with CIB microbiome data###################
path <- setwd("C:/Users/arial/Desktop/Ch.3 Microbiome & popgen/CIB-ubiome-popgen-code/2016-2018_fastq/")
list.files(path)

#Forward and Reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

###Inspect read quality profiles
#Forward reads
plotQualityProfile(fnFs)

#Reverse reads
plotQualityProfile(fnRs)

###Filter and trim
#Assign filenames for the filtered fastq.gz files and place filtered files in filtered/sub-directory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. 
#maxEE parameter sets the max number of “expected errors” allowed in a read, which is a better than simply averaging quality scores.
#AVC did no trim left, used tuncQ=10 to truncate quality scores above 10
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncQ=10, maxN=0, maxEE=c(2,2), rm.phix=TRUE, compress=TRUE, multithread=FALSE) 
head(out)

###Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#Visualize estimated error rates
plotErrors(errF, nominalQ=TRUE)

###Sample inference
#apply the core sample inference algorithm to the filtered and trimmed sequence data
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Inspect the first sample
dadaFs[[1]]

###Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#Inspect the merger data.frame from the first sample
head(mergers[[1]])

###Construct ASV sequence table
seqtab <- makeSequenceTable(mergers)
#samples and ASVs
dim(seqtab) 

#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

###Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#view results
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#save results
write.csv(seqtab.nochim, "C:/Users/arial/Desktop/Ch.3 Microbiome & popgen/CIB-ubiome-popgen-code/seqtab.nochim.csv")


###Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#save results
write.csv(track, "C:/Users/arial/Desktop/Ch.3 Microbiome & popgen/CIB-ubiome-popgen-code/track.csv")


###Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/arial/Desktop/Ch.3 Microbiome & popgen/CIB-ubiome-popgen-code/silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#save taxa results
write.csv(taxa, "C:/Users/arial/Desktop/Ch.3 Microbiome & popgen/CIB-ubiome-popgen-code/taxa.csv")
################################################################################


###DADA2 Handoff to Phyloseq
#Set path
path <- setwd("C:/Users/arial/Desktop/Ch.3 Microbiome & popgen/CIB-ubiome-popgen/Data/")
list.files(path)

#Read in saved dada2 files
asv_mat <- as.matrix(read.csv("seqtab.nochim.csv", row.names=1))
taxa_mat <- as.matrix(read.csv("taxa.csv", row.names=1))
samples <- read.csv("VanCise_ubiome_metadata.csv", row.names=1)

#Create phyloseq components
asv <- otu_table(asv_mat, taxa_are_rows=FALSE)
taxa <- tax_table(taxa_mat)
samples <- sample_data(samples)
samples$year <- as.factor(samples$year)

#Create phyloseq object
physeq <- phyloseq(asv, taxa, samples)
ntaxa(physeq)

#examine sample taxa
plot_bar(physeq, fill = "Kingdom")

#Remove eukaryotic sequences
physeq <- subset_taxa(physeq, Kingdom != "Eukaryota")
ntaxa(physeq)
plot_bar(physeq, fill = "Kingdom")

#remove samples with <10,000 reads
physeq <- prune_samples(sample_sums(physeq) >= 10000, physeq)
ntaxa(physeq)
#no samples had <10,000 reads

#create table to see number per phyla
table(tax_table(physeq)[, "Phylum"], exclude = NULL)

#Explore alpha diversity among samples
plot_richness(physeq, measures=c("Shannon", "Simpson"))

#Explore alpha diversity among variables
plot_richness(physeq, x="AgeClass", measures=c("Shannon", "Simpson"))
plot_richness(physeq, x="sex", measures=c("Shannon", "Simpson"))
plot_richness(physeq, x="year", measures=c("Shannon", "Simpson"))
plot_richness(physeq, x="Location", measures=c("Shannon", "Simpson"))

#Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(physeq, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color="AgeClass", title="Bray NMDS")
plot_ordination(ps.prop, ord.nmds.bray, color="sex", title="Bray NMDS")
plot_ordination(ps.prop, ord.nmds.bray, color="year", title="Bray NMDS")
plot_ordination(ps.prop, ord.nmds.bray, color="Location", title="Bray NMDS")
#possible clusters with year and location


###Look at relative abundances
#at the kingdom level
ps_kingdom <- tax_glom(physeq, "Kingdom", NArm = TRUE)

#Transform Taxa counts to relative abundance
ps_kingdom_relabun <- transform_sample_counts(ps_kingdom, function(OTU) OTU/sum(OTU) * 100)
taxa_abundance_table_kingdom <- psmelt(ps_kingdom_relabun)

#plot
ggplot(data=taxa_abundance_table_kingdom, aes(x =Sample, y = Abundance, fill = Kingdom)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Relative Abundance") +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 12)) +
  scale_y_continuous(expand=c(0,0))


#at the phylum level
ps_phylum <- tax_glom(physeq, "Phylum", NArm = TRUE)

#Transform Taxa counts to relative abundance
ps_phylum_relabun <- transform_sample_counts(ps_phylum, function(OTU) OTU/sum(OTU) * 100)
taxa_abundance_table_phylum <- psmelt(ps_phylum_relabun)

#plot
ggplot(data=taxa_abundance_table_phylum, aes(x =Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Relative Abundance") +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
       axis.text.y = element_text(size = 12),
       legend.text = element_text(size = 10),
       strip.text = element_text(size = 12)) +
  scale_y_continuous(expand=c(0,0)) +
  facet_grid(~ Location, scales = "free") 

#explore boxplots
taxa_abundance_table_phylum %>% 
  ggplot(aes(x =Phylum, y = Abundance, fill = Phylum)) +
  geom_boxplot() +
  labs(x = "",
       y = "Relative Abundance",
       title = "Phylum Relative Abundance") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12)) +
  scale_y_continuous(expand=c(0,0)) +
  facet_grid(~ Location, scales = "free") 



#It is more convenient to use short names for our ASVs (e.g. ASV21) rather than the full DNA sequence when working with
#some of the tables and visualizations from phyloseq, but we want to keep the full DNA sequences for other purposes like
#merging with other datasets or indexing into reference databases. For that reason we’ll store the DNA sequences of our
#ASVs in the refseq slot of the phyloseq object, and then rename our taxa to a short string. That way, the short new taxa
#names will appear in tables and plots, and we can still recover the DNA sequences corresponding to each ASV as needed with refseq(ps).
# dna <- Biostrings::DNAStringSet(taxa_names(physeq))
# names(dna) <- taxa_names(physeq)
# physeq <- merge_phyloseq(physeq, dna)
# taxa_names(physeq) <- paste0("ASV", seq(ntaxa(physeq)))
# physeq

#Calculate Bray-Curtis distance
bray_dist <- distance(physeq, method = "bray")

#save as matrix
bray_dist_mat <- as.matrix(bray_dist)

#subtract dissimilarity values by 1 to get similarity
bray_sim_mat <- 1-bray_dist_mat

#remove exponential notation to get fixed decimal format
bray_sim_mat <- format(bray_sim_mat, scientific = FALSE)


#save microbial similarity results
write.csv(bray_sim_mat, "C:/Users/arial/Desktop/Ch.3 Microbiome & popgen/CIB-ubiome-popgen/ubiome_matrix.csv")












