#This script scrapes the three publicly available datasets off of GEO and prepares them for module preservation analysis
#These datasets are cleaned such that only genes found in our original data will be included in them 

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
colors <- vector(mode = "list", length = 4)
samples <- vector(mode = "list", length = 4)
genes <- vector(mode = "list", length = 4)
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

dcolon.1 <- dcolon[, genes[[1]]]
dcolon.2 <- dcolon[, genes[[2]]]
dcolon.3 <- dcolon[, genes[[4]]]

pcolon.1 <- pcolon[, genes[[1]]]
pcolon.2 <- pcolon[, genes[[3]]]
pcolon.3 <- pcolon[, genes[[4]]]

#Distal colon vs GSE131032----
setLabels <- c("Distal", "GSE131032")
multiExpr <- list(Distal = list(data = dcolon.1), GSE131032 = list(data = dfs[[1]]))
multiColor <- list(Distal = colors[[1]], GSE131032 = colors[[1]])
nSets <- 2

system.time({
  mp.d1 = modulePreservation(multiExpr, multiColor,
                             referenceNetworks = 1,
                             nPermutations = 200,
                             randomSeed = 1,
                             verbose = 5,
                             networkType = "signed hybrid",
                             parallelCalculation = F)
})
#Distal colon vs GSE168053 (distal)----
setLabels <- c("Distal", "GSE168053_dist")
multiExpr <- list(Distal = list(data = dcolon.2), GSE131032 = list(data = dfs[[2]]))
multiColor <- list(Distal = colors[[2]], GSE131032 = colors[[2]])
nSets <- 2

system.time({
  mp.d2 = modulePreservation(multiExpr, multiColor,
                             referenceNetworks = 1,
                             nPermutations = 200,
                             randomSeed = 1,
                             verbose = 5,
                             networkType = "signed hybrid",
                             parallelCalculation = F)
})
#Distal colon vs GSE210405----
setLabels <- c("Distal", "GSE210405")
multiExpr <- list(Distal = list(data = dcolon.3), GSE131032 = list(data = dfs[[4]]))
multiColor <- list(Distal = colors[[4]], GSE131032 = colors[[4]])
nSets <- 2

system.time({
  mp.d3 = modulePreservation(multiExpr, multiColor,
                             referenceNetworks = 1,
                             nPermutations = 200,
                             randomSeed = 1,
                             verbose = 5,
                             networkType = "signed hybrid",
                             parallelCalculation = F)
})
#Proximal colon vs GSE131032----
setLabels <- c("Proximal", "GSE131032")
multiExpr <- list(Proximal = list(data = pcolon.1), GSE131032 = list(data = dfs[[1]]))
multiColor <- list(Proximal = colors[[1]], GSE131032 = colors[[1]])
nSets <- 2

system.time({
  mp.p1 = modulePreservation(multiExpr, multiColor,
                             referenceNetworks = 1,
                             nPermutations = 200,
                             randomSeed = 1,
                             verbose = 5,
                             networkType = "signed hybrid",
                             parallelCalculation = F)
})
#Proximal colon vs GSE168053 (proximal)----
setLabels <- c("Proximal", "GSE168053_dist")
multiExpr <- list(Proximal = list(data = pcolon.2), GSE131032 = list(data = dfs[[3]]))
multiColor <- list(Proximal = colors[[3]], GSE131032 = colors[[3]])
nSets <- 2

system.time({
  mp.p2 = modulePreservation(multiExpr, multiColor,
                             referenceNetworks = 1,
                             nPermutations = 200,
                             randomSeed = 1,
                             verbose = 5,
                             networkType = "signed hybrid",
                             parallelCalculation = F)
})
#Proximal colon vs GSE210405----
setLabels <- c("Proximal", "GSE210405")
multiExpr <- list(Proximal = list(data = pcolon.3), GSE131032 = list(data = dfs[[4]]))
multiColor <- list(Proximal = colors[[4]], GSE131032 = colors[[4]])
nSets <- 2

system.time({
  mp.p3 = modulePreservation(multiExpr, multiColor,
                             referenceNetworks = 1,
                             nPermutations = 200,
                             randomSeed = 1,
                             verbose = 5,
                             networkType = "signed hybrid",
                             parallelCalculation = F)
})
