#This script scrapes the three publicly available datasets off of GEO and prepares them for module preservation analysis
#These datasets are cleaned such that only genes found in our original data will be included in them 
#These cleaned expression matrices are saved for the preservation analysis and later calculation of eigengenes, etc

library(openxlsx)
library(tidyverse)
library(WGCNA)
library(plyr)

enableWGCNAThreads()

load("colon_RNAseq/data_objects/consensus_network.RData")
load("colon_RNAseq/consensus_WGCNA_input.RData")

#Scrape the first dataset from GEO
gse1 <- gzcon(url("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE131032&format=file&file=GSE131032%5Fkallisto%5Fcounts%2Ecsv%2Egz"))
gse1 <- readLines(gse1)
gse1 <- read.csv(textConnection(gse1))

#Get second dataset 
gse2 <- gzcon(url("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE168053&format=file&file=GSE168053%5Fabundances%2Ecsv%2Egz"))
gse2 <- readLines(gse2)
gse2 <- read.csv(textConnection(gse2))

#Get third dataset
gse3 <- read.xlsx("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE210405&format=file&file=GSE210405%5Fnormalized%5Fcount%5Fmatrix%2Exlsx")

genes <- colnames(multiExpr[[1]]$data)

rownames(gse1) <- NULL
rownames(gse2) <- NULL
rownames(gse3) <- NULL

gse1 <- gse1 %>%
  column_to_rownames(var = "X") %>%
  t()

gse2_distal <- gse2 %>%
  na.omit() %>%
  plyr::ddply("ext_gene", numcolwise(sum)) %>%
  dplyr::select(-entrez_id) %>%
  column_to_rownames(var = "ext_gene") %>%
  dplyr::select(matches("dist")) %>%
  t()

gse2_proximal <- gse2 %>%
  na.omit() %>%
  plyr::ddply("ext_gene", numcolwise(sum)) %>%
  dplyr::select(-entrez_id) %>%
  column_to_rownames(var = "ext_gene") %>%
  dplyr::select(matches("prox")) %>%
  t()

gse3 <- gse3 %>%
  column_to_rownames(var = "X1") %>%
  t()

dfs <- list(gse1, gse2_distal, gse2_proximal, gse3)
names(dfs) <- c("gse1", "gse2_distal", "gse2_proximal", "gse3")

for (i in seq_along(dfs)){
  dfs[[i]] <- as.data.frame(dfs[[i]]) %>%
    dplyr::select(intersect(genes, colnames(dfs[[i]])))
  gsg <- goodSamplesGenes(dfs[[i]])
  dfs[[i]] <- dfs[[i]][gsg$goodSamples, gsg$goodGenes]
}

key <- data.frame(color = moduleColors,
                  gene = genes)
keys <- vector(mode = "list", length = 4)
names(keys) <- names(dfs)
colors <- vector(mode = "list", length = 4)
names(colors) <- names(dfs)
samples <- vector(mode = "list", length = 4)
names(samples) <- names(dfs)
genes <- vector(mode = "list", length = 4)
names(genes) <- names(dfs)
for (i in seq_along(keys)){
  #Make named color vector, which identifies module members
  keys[[i]] <- key[key$gene %in% colnames(dfs[[i]]), ]
  colors[[i]] <- keys[[i]]$color
  names(colors[[i]]) <- keys[[i]]$gene
  
  #Pull out genes and samples
  genes[[i]] <- colnames(dfs[[i]])
  samples[[i]] <- rownames(dfs[[i]])
  
  #Make sure public data sets are numeric expresion matrices
  dfs[[i]] <- apply(dfs[[i]], 2, as.numeric)
  rownames(dfs[[i]]) <- samples[[i]]
}

dcolon <- multiExpr[[1]]$data
pcolon <- multiExpr[[2]]$data

dcolon.1 <- dcolon[, genes$gse1]
dcolon.2 <- dcolon[, genes$gse2_distal]
dcolon.3 <- dcolon[, genes$gse3]

pcolon.1 <- pcolon[, genes$gse1]
pcolon.2 <- pcolon[, genes$gse2_distal]
pcolon.3 <- pcolon[, genes$gse3]

save(dcolon.1, dcolon.2, dcolon.3,
     pcolon.1, pcolon.2, pcolon.3,
     colors, keys, samples, genes, 
     dfs,
     file = "colon_RNAseq/data_objects/module_preservation_input.RData")