rm(list=ls()) #clean environment
#dada2_singleread
library(dada2); packageVersion("dada2")
#setup paths and folders
setwd("~/Documents/Microbiome_data/PRJNA589343_trimmed/")
path <- "~/Documents/Microbiome_data/PRJNA589343_trimmed/"
filtpath <- file.path(path, "filtered")
dir.create(filtpath)
fns <- list.files(path, pattern="_trim.fastq.gz")
filts <- file.path(filtpath, fns)
#check quality of reads
plotQualityProfile(fns[1:4])
#filter and trim each sample file
filterAndTrim(fns, filts, 
              truncLen=145, maxEE=2, truncQ=2, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=FALSE)
#check quality of reads
plotQualityProfile(filts[1:4])
#Sample Inference

#file parsing
filtpath <- "filtered/" 
fns <- list.files(filtpath, pattern="_trim.fastq.gz", full.names = TRUE) # CHANGE if different file extensions
sample.names <- sapply(strsplit(basename(filts), "_"), `[`, 1)
names(filts) <- sample.names
#Learn error rates
set.seed(100)
err <- learnErrors(filts, nreads = 2e6, multithread=TRUE, randomize=TRUE)
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
write.csv(nonchim, "Bimera_output.csv")

#assign taxonomy
taxa_silva <- assignTaxonomy(seqtab.nochim, "~/Documents/Microbiome_data/Michael-and-Ben-Repo/Datasets/silva_nr_v132_train_set.fa.gz", tryRC=TRUE, multithread=TRUE)
#unname(head(taxa))
write.csv(taxa_silva, "taxa_PRJNA589343.csv")
saveRDS(taxa_silva, "taxa_PRJNA589343.rds")
#unname(head(taxa_gg))
unname(head(taxa_silva))

#write to disk
saveRDS(seqtab.nochim, "seqtab_PRJNA589343.rds")
write.csv(seqtab.nochim, "seqtab_PRJNA589343.csv")
seqtab.nochim.trans = t(seqtab.nochim)
write.csv(seqtab.nochim.trans, "seqtab_transpose_PRJNA589343.csv")

# Assign taxonomy and write to disk

