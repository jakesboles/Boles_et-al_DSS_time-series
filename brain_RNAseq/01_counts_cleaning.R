library(DESeq2)
library(tidyverse)
library(gtools)
library(edgeR)

#Grab raw counts file from GEO's https link and convert to data frame----

con <- gzcon(url("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE239820&format=file&file=GSE239820%5Fraw%5Fcounts%2Ecsv%2Egz"))
csv <- readLines(con)
counts <- read.csv(textConnection(csv))

counts <- counts[, mixedsort(colnames(counts))] %>%
  column_to_rownames(var = "X") %>%
  mutate_all(function(x) as.numeric(x)) %>%
  dplyr::select(JB1:JB265) #these are the brain samples 

#Using CPM from the edgeR package to filter genes with low counts----
cpm <- cpm(counts)

col1sum <- sum(counts[, 1]) / 1000000 #these two lines ensure CPM was calculated correctly
counts[1, 1]/col1sum                  #this line should return 0  

threshold <- cpm > 0.5
keep <- rowSums(threshold) >= 8
summary(keep)
counts <- counts[keep, ]
dim(counts) #should be left with 21603 genes (rows) and 265 samples (columns)

#Load sample data, create DESeq2 objects and save variance-stabilized data----
brain_key <- read.csv("brain_RNAseq/RNAseq_brain_metadata.csv")
brain_key <- brain_key %>%
  mutate(tissue = factor(tissue),
         group = factor(group,
                        levels = c("Untreated", "5d_DSS", "7d_DSS", "7d_DSS_2d_H2O",
                                   "7d_DSS_5d_H2O", "7d_DSS_7d_H2O",
                                   "7d_DSS_14d_H2O")))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = brain_key,
                              design = ~ tissue + group)
dds

dds <- DESeq(dds)

write.csv(assay(vsd),                            #write this to a CSV file for next step
          file = "brain_RNAseq/vst_counts_1.csv")