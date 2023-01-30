rm(list=ls()) #clean environment
#dada2_singleread
library(dada2); packageVersion("dada2")
#setup paths and folders

#setwd("~/Desktop/Disk1/Saad/190715_JDS/")
path <- "/home/agronah/Documents/Microbiome_data/PRJNA168470"
filtpath <- file.path(path, "filtered")
dir.create(filtpath)
fns <- list.files(path, pattern=".fastq.gz")
filts <- file.path(filtpath, fns)
#check quality of reads
#plotQualityProfile(fns[1:4])
#filter and trim each sample file  truncLen=145,
#truncLen=200, maxEE=5,
out <- filterAndTrim(fns, filts, 
                     truncLen=200, maxEE=5, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=FALSE)

#head(out) 
#check quality of reads
plotQualityProfile(filts[1:4])
#Sample Inference

#file parsing
filtpath <- "filtered/" 
fns <- list.files(filtpath, pattern=".fastq.gz", full.names = TRUE) # CHANGE if different file extensions
sample.names <- sapply(strsplit(basename(filts), "_"), `[`, 1)
names(filts) <- sample.names
#Learn error rates
set.seed(100)
err <- learnErrors(filts, nreads = NULL, multithread=TRUE, randomize=TRUE)

plotErrors(err, nominalQ=TRUE)

# Infer sequence variants
dds <- vector("list", length(sample.names))
names(dds) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derep <- derepFastq(filts[[sam]])
  dds[[sam]] <- dada(derep, err=err, multithread=TRUE)
}
#construct sequence table
seqtab <- makeSequenceTable(dds)
#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
nonchim = sum(seqtab.nochim)/sum(seqtab)
write.csv(nonchim, "Bimera_output_PRJNA168470.csv")

#assign taxonomy
taxa_silva <- assignTaxonomy(seqtab.nochim, "~/Documents/Microbiome_data/Michael-and-Ben-Repo/Datasets/silva_nr_v132_train_set.fa.gz", tryRC=TRUE, multithread=TRUE)

write.csv(taxa_silva, "taxa_PRJNA168470.csv")
saveRDS(taxa_silva, "taxa_PRJNA168470.rds")
#unname(head(taxa_gg))
unname(head(taxa_silva))

#write to disk
saveRDS(seqtab.nochim, "seqtab_R2_v2v3_PRJNA168470.rds")
write.csv(seqtab.nochim, "seqtab_R2_v2v3_PRJNA168470.csv")
seqtab.nochim.trans = t(seqtab.nochim)
write.csv(seqtab.nochim.trans, "seqtab_trans_v2v3_PRJNA168470.csv")

# Assign taxonomy and write to disk

