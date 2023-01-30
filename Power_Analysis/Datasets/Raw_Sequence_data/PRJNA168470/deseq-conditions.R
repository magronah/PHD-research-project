Data_PRJNA168470 = read.csv("~/Documents/Microbiome_data/Michael-and-Ben-Repo/Datasets/PRJNA168470/seqtab_nochim_lessHOST_PRJNA168470.csv", row.names = 1) 
names <- colnames(Data_PRJNA168470)
groups <- read.csv("groups_PRJNA168470.csv")
grps <- as.factor(groups$x)
condition <- data.frame(names, grps)

met_data <- read.csv("metadata_PRJNA168470.txt")

g<-substring(met_data$Sample.Name, 1,1)
stopifnot(all(met_data$g==file$groups))
stopifnot(all(met_data$Run==names))

write.csv(condition, file="Deseq_condition_PRJNA168470.csv")
conditions_PRJNA168470 <- read.csv("~/Documents/Microbiome_data/Michael-and-Ben-Repo/Datasets/PRJNA168470/Deseq_condition_PRJNA168470.csv", header= T, row.names = 1) 
str(conditions_PRJNA168470)



