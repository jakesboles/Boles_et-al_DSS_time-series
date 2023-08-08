library(DESeq2)
library(tidyverse)
library(gtools)
library(edgeR)
library(WGCNA)

#Grab raw counts file from GEO's https link and convert to data frame

con <- gzcon(url("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE239820&format=file&file=GSE239820%5Fraw%5Fcounts%2Ecsv%2Egz"))
csv <- readLines(con)
counts <- read.csv(textConnection(csv))

counts <- counts[, mixedsort(colnames(counts))] %>%
  column_to_rownames(var = "X") %>%
  mutate_all(function(x) as.numeric(x)) %>%
  dplyr::select(JB266:JB395) #these are the colon samples 

cpm <- cpm(counts)

col1sum <- sum(counts[, 1]) / 1000000 #these two lines ensure CPM was calculated correctly
counts[1, 1]/col1sum                  #this line should return 0  

threshold <- cpm > 0.5
keep <- rowSums(threshold) >= 8
summary(keep)
counts <- counts[keep, ]
dim(counts) #should be left with 20604 genes (rows) and 130 samples (columns)
