rm(list=ls())

setwd("~/Documents/Microbiome_data/PRJNA589343_trimmed")
path <- "~/Documents/Microbiome_data/PRJNA589343_trimmed"

# load libraries for ordination and clustering
library("cluster")
library(data.table)
library("phyloseq")
packageDescription("phyloseq")$Version
library("ggplot2")
library("plyr")
library("grid")
library("ape")
library("phangorn")
library("phytools")
library("vegan")
 
# Read in OTU table
otu_df = read.csv("seqtab_transpose_PRJNA589343.csv", row.names = 1) 
seqs = rownames(otu_df)
rownames(otu_df) = NULL
 
#Read in the taxfile 
tax_df = read.csv("taxa_PRJNA589343.csv", row.names = 1)
#View(tax_df)

#Be sure that the row names of taxa and OTU table  matches
stopifnot(all(seqs == rownames(tax_df)))

rownames(tax_df) = NULL #remove row name for taxa  

#Read in the Metadata
met_df = read.csv("metadata_PRJNA589343.txt") 
View(met_df) 
#dim(met_df)

#Extracting group information
gp <- met_df$Sample.Name
groups <- substring(gp,1,3)
unique(groups)
sum(groups == "HCs")
sum(groups == "ASD")

kk <- data.frame(met_df$Run, substring(met_df$Sample.Name, 1,3))

colnames(kk) <- c("Samples", "Groups")
write.csv(kk,file= "groups_names_PRJNA589343.csv", row.names= FALSE)
#write.csv(groups,file= "groups_PRJNA589343.csv", row.names= FALSE)

#Some extra cleaning
met_df$Run = gsub('-', '.', as.character(met_df$Run))
met_df$Run = gsub(' ', '.', met_df$Run)
met_df$Run = as.factor(met_df$Run)
rownames(met_df) = met_df$Run   # make sure the rownames in the metadata data frame are the same as the row (or column, whichever is samples) in the otu data frame

#create the phyloseq object
dat = phyloseq(otu_table(otu_df, taxa_are_rows = TRUE), # or FALSE if false
               tax_table(as.matrix(tax_df)),
               sample_data(met_df))

#preprocessing 
#write ASV sequences to file
write.csv(seqs, file="ASV__PRJNA589343_sequences.csv")
 
#try writing fasta file with new package **this works!
#seqs.fasta = dataframe2fas(seqs, file="ASVseqs.fasta")

#total reads per sample and distribution
sdt = data.table(as(sample_data(dat), "data.frame"),
                 TotalReads = sample_sums(dat), keep.rownames = TRUE)
#View(sdt) 

setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
pSeqDepth + facet_wrap(~Sample)

#simple sums and write to file
samplesums = sort(sample_sums(dat))
write.csv(samplesums, file="readcounts_per_sample_sorted_PRJNA589343.csv")
plot(samplesums)

#filter table of non-bacterial taxa
dat_lessHOST = subset_taxa(dat, Kingdom=="Bacteria", Family!="Mitochondria")
write.csv(otu_table(dat_lessHOST),file='seqtab_nochim_lessHOST__PRJNA589343.csv')
write.csv(tax_table(dat_lessHOST),file='taxa_JJML_lessHOST__PRJNA589343.csv')

#remove less than 0.001% ?
rel_abun_all = transform_sample_counts(dat, function(x) x/sum(x))
rel_abun_all_prune = prune_taxa(taxa_sums(rel_abun_all) > 0.001,
                                rel_abun_all)

write.csv(otu_table(rel_abun_all_prune),file='relative_abundance__PRJNA589343.csv')
rel_abun_PRJNA589343 <- read.csv("relative_abundance__PRJNA589343.csv")
#View(rel_abun_PRJNA589343)

#remove less than 5 reads
dat_lessHOST_n5 = prune_taxa(taxa_sums(dat_lessHOST)>5, dat_lessHOST)

#rarefy
#min_lib <- min(sample_sums(dat_lessHOST)) 
dat_r =rarefy_even_depth(dat_lessHOST,sample.size=5000,verbose=FALSE,replace = FALSE)
dat_r1000=rarefy_even_depth(dat_lessHOST, sample.size=1000, verbose=FALSE, replace=FALSE)
write.csv(otu_table(dat_r),file="rarefy_otu_PRJNA589343.csv")
