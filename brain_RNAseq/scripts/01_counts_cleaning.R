#In this file, we scrape the raw counts file from our GEO accession page (linked in the README)
#We first filter on counts-per-million using the edgeR package, then we use DESeq2 and WGCNA to remove bad samples and genes
#Samples are filtered based on distance in a dendrogram, and genes are further excluded on the basis of having zero variance or too many zero counts
#After each gene/sample removal stage, we adjust the raw counts to reflect the cleaned data and re-normalize the data using DESeq2 to produce variance-stabilized counts.
#At the end, we will have a .csv file containing the variance-stabilized counts, which will be the input to the WGCNA, and a DESeq2 object, which we can use for differential expression analyses.

library(DESeq2)
library(tidyverse)
library(gtools)
library(edgeR)
library(WGCNA)

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
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = T)
vsd <- assay(vsd)

#Using WGCNA's system to discard low-quality samples/genes----
#Create list of expression matrices reflecting the four "networks" 
ctx <- as.data.frame(vsd) %>%
  dplyr::select(JB1:JB69)
mb <- as.data.frame(vsd) %>%
  dplyr::select(JB70:JB133)
hpc <- as.data.frame(vsd) %>%
  dplyr::select(JB134:JB198)
str <- as.data.frame(vsd) %>%
  dplyr::select(JB199:JB265)

nSets <- 4
setLabels <- c("Cortex", "Midbrain", "Hippocampus", "Striatum")
shortLabels <- c("CTX", "MB", "HPC", "STR")

multiExpr <- vector(mode = "list", length = nSets)
multiExpr[[1]] <- list(data = as.data.frame(t(ctx)))
rownames(multiExpr[[1]]$data) <- names(ctx)

multiExpr[[2]] <- list(data = as.data.frame(t(mb)))
rownames(multiExpr[[2]]$data) <- names(mb)

multiExpr[[3]] <- list(data = as.data.frame(t(hpc)))
rownames(multiExpr[[3]]$data) <- names(hpc)

multiExpr[[4]] <- list(data = as.data.frame(t(str)))
rownames(multiExpr[[4]]$data) <- names(str)

exprSize <- checkSets(multiExpr)
exprSize

#Mark and remove low quality genes
gsg <- goodSamplesGenesMS(multiExpr, verbose = 10)
gsg$allOK

if (!gsg$allOK){
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}
exprSize

#Mark and remove low quality samples with the assistance of a dendrogram 
sampleTrees <- list()
for (set in 1:nSets){
  sampleTrees[[set]] <- hclust(dist(multiExpr[[set]]$data), method = "average")
}
pdf(file = "brain_RNAseq/plots/sample_clustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets){
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7)
}
dev.off()

cutHeights <- c(150, 140, 100, 140) #cut height to separate "outliers" from the rest
pdf(file = "brain_RNAseq/plots/sample_clustering_with_cuts.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets){
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
  abline(h=cutHeights[set], col = "red");
}
dev.off()

for (set in 1:nSets){
  # Find clusters cut by the line
  labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
  # Keep the largest one (labeled by the number 1)
  keep = (labels==1)
  multiExpr[[set]]$data = multiExpr[[set]]$data[keep, ]
}

collectGarbage()
# Check the size of the leftover data
exprSize <- checkSets(multiExpr)
exprSize

#Remove low quality genes/samples from raw data to then re-normalize, etc.
gene_index <- rownames(counts) %in% colnames(multiExpr[[1]]$data)
counts <- counts[gene_index, ]

ctx_counts <- counts %>% #separate into regions to ease sample removal
  dplyr::select(JB1:JB69)
mb_counts <- counts %>%
  dplyr::select(JB70:JB133)
hpc_counts <- counts %>%
  dplyr::select(JB134:JB198)
str_counts <- counts %>%
  dplyr::select(JB199:JB265)

ctx_sample_index <- colnames(ctx_counts) %in% rownames(multiExpr[[1]]$data) 
mb_sample_index <- colnames(mb_counts) %in% rownames(multiExpr[[2]]$data)
hpc_sample_index <- colnames(hpc_counts) %in% rownames(multiExpr[[3]]$data)
str_sample_index <- colnames(str_counts) %in% rownames(multiExpr[[4]]$data)

ctx_counts <- ctx_counts[, ctx_sample_index]
mb_counts <- mb_counts[, mb_sample_index]
hpc_counts <- hpc_counts[, hpc_sample_index]
str_counts <- str_counts[, str_sample_index]

counts <- cbind(ctx_counts, mb_counts, hpc_counts, str_counts) #rebuild raw counts

#Update metadata frame 
brain_key <- brain_key[brain_key$FileID %in% colnames(counts), ]

#Re-normalize and variance stabilize with DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = brain_key,
                              design = ~ tissue + group)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = T)
vsd <- assay(vsd)

#Recheck gene quality with the same steps as above
ctx <- vsd[, colnames(vsd) %in% brain_key[grep("cortex", brain_key$tissue), ]$FileID]
mb <- vsd[, colnames(vsd) %in% brain_key[grep("midbrain", brain_key$tissue), ]$FileID]
hpc <- vsd[, colnames(vsd) %in% brain_key[grep("hippocampus", brain_key$tissue), ]$FileID]
str <- vsd[, colnames(vsd) %in% brain_key[grep("striatum", brain_key$tissue), ]$FileID]  

multiExpr <- vector(mode = "list", length = nSets)
multiExpr[[1]] <- list(data = as.data.frame(t(ctx)))
rownames(multiExpr[[1]]$data) <- names(ctx)

multiExpr[[2]] <- list(data = as.data.frame(t(mb)))
rownames(multiExpr[[2]]$data) <- names(mb)

multiExpr[[3]] <- list(data = as.data.frame(t(hpc)))
rownames(multiExpr[[3]]$data) <- names(hpc)

multiExpr[[4]] <- list(data = as.data.frame(t(str)))
rownames(multiExpr[[4]]$data) <- names(str)

exprSize <- checkSets(multiExpr)
exprSize

gsg <- goodSamplesGenesMS(multiExpr, verbose = 10)
gsg$allOK #this flags 4 more genes for removal 

if (!gsg$allOK){
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}
exprSize

#Remove flagged genes from raw data (again)
gene_index <- rownames(counts) %in% colnames(multiExpr[[1]]$data)
counts <- counts[gene_index, ]

#Normalize (again)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = brain_key,
                              design = ~ tissue + group)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = T)
vsd <- assay(vsd)

#Check for bad genes (again)
ctx <- vsd[, colnames(vsd) %in% brain_key[grep("cortex", brain_key$tissue), ]$FileID]
mb <- vsd[, colnames(vsd) %in% brain_key[grep("midbrain", brain_key$tissue), ]$FileID]
hpc <- vsd[, colnames(vsd) %in% brain_key[grep("hippocampus", brain_key$tissue), ]$FileID]
str <- vsd[, colnames(vsd) %in% brain_key[grep("striatum", brain_key$tissue), ]$FileID]  

multiExpr <- vector(mode = "list", length = nSets)
multiExpr[[1]] <- list(data = as.data.frame(t(ctx)))
rownames(multiExpr[[1]]$data) <- names(ctx)

multiExpr[[2]] <- list(data = as.data.frame(t(mb)))
rownames(multiExpr[[2]]$data) <- names(mb)

multiExpr[[3]] <- list(data = as.data.frame(t(hpc)))
rownames(multiExpr[[3]]$data) <- names(hpc)

multiExpr[[4]] <- list(data = as.data.frame(t(str)))
rownames(multiExpr[[4]]$data) <- names(str)

exprSize <- checkSets(multiExpr)
exprSize

gsg <- goodSamplesGenesMS(multiExpr, verbose = 10)
gsg$allOK #all good now! 

#Scrape uploaded VST-counts from our GEO submission to compare
con <- gzcon(url("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE239820&format=file&file=GSE239820%5Fbrain%5Fclean%5FVST%5Fcounts%2Ecsv%2Egz"))
lines <- readLines(con)
vsd0 <- read.csv(textConnection(lines))
vsd0 <- vsd0 %>%
  column_to_rownames(var = "X") %>%
  as.matrix()

all.equal(vsd, vsd0)
#This should return TRUE, indicating the matrix we've created here is identical to that used in the paper 

#Save cleaned VST counts for later use
write.csv(vsd, file = "brain_RNAseq/csv_outputs/vst_counts.csv",
          col.names = T, row.names = T)
save(dds, brain_key, vsd,
     file = "brain_RNAseq/data_objects/deseq2_objects.RData")

