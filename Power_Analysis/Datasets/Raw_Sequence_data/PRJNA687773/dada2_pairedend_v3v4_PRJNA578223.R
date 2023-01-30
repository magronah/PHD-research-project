rm(list=ls())
library(dada2); packageVersion("dada2")
path = "~/Documents/Microbiome_data/PRJNA578223"  

list.files(path)
# Sort ensures forward/reverse reads are in same order

fnFs <- sort(list.files(path, pattern="1.fastq.gz"))
fnRs <- sort(list.files(path, pattern="2.fastq.gz"))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

#check quality profile
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])

#filter and trim
#define filenames
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

#filter forward and reverse reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,220),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, verbose=TRUE)  
head(out)

#check quality of reads
png(file="QualityProfile_ForwardRead",
    width=800, height=450)
plotQualityProfile(filtFs[1:4])
dev.off()

png(file="QualityProfile_ReverseRead",
    width=800, height=450)
plotQualityProfile(filtRs[1:4])
dev.off()
#learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#visualize estimated error rates
png(file="ErrorPlot_ForwardRead",
    width=800, height=450)
plotErrors(errF, nominalQ=TRUE)
dev.off()

png(file="ErrorPlot_ReverseRead",
    width=800, height=450)
plotErrors(errR, nominalQ=TRUE)
dev.off()
#dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#infer sequence variants
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#inspect dada-class object
#dadaFs[[1]]
#dadaRs[[1]]

#merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE, minOverlap=5)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#make sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) 

#SRR13306642 
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
nonchim = sum(seqtab.nochim)/sum(seqtab)

write.csv(nonchim, "Bimera_output_PRJNA578223.csv")

#track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "read_track_PRJNA578223.csv")

#assign taxonomy
taxa_silva <- assignTaxonomy(seqtab.nochim, "~/Documents/Microbiome_data/Michael-and-Ben-Repo/Datasets/silva_nr_v132_train_set.fa.gz", tryRC=TRUE, multithread=TRUE)
#unname(head(taxa))

write.csv(taxa_silva, "taxa_silva132_v3v4_PRJNA578223.csv")
saveRDS(taxa_silva, "taxa_silva132_v3v4_PRJNA578223.rds")

#unname(head(taxa_gg))
unname(head(taxa_silva))

#write to disk
write.csv(seqtab.nochim, "seqtab_nochim_v3v4_PRJNA578223.csv")
saveRDS(seqtab.nochim, "seqtab_nochim_v3v4_PRJNA578223.rds")
seqtab.nochim.trans = t(seqtab.nochim)
write.csv(seqtab.nochim.trans, "seqtab_nochim_transposed_v3v4_PRJNA578223.csv")

#file clean-up
#system("rm filtered/*.fastq.gz")
#system("rm *_trim*.fastq.gz*")
