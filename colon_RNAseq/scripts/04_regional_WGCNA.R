#This script creates weighted co-expression networks from both colon segments individually, using parameters identified in 03
#The only purpose of this analysis is to later evaluate how similar the consensus modules are to segment-specific modules

library(WGCNA)

load("colon_RNAseq/consensus_WGCNA_input.RData")

plots <- "colon_RNAseq/plots/"
data <- "colon_RNAseq/data_objects/"

options(stringsAsFactors = F)
enableWGCNAThreads()

message("Starting distal colon")
#Distal colon----
datExpr <- multiExpr[[1]]$data

adjacency <- adjacency(datExpr, power = softPower, type = "signed hybrid")

TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")

sizeGrWindow(12,9)
pdf(paste0(plots, "distal_colon_gene_clustering.pdf", height = 9, width = 12))
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity\nDistal colon",
     labels = F, hang = 0.04)
dev.off()

minModuleSize <- 30

dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = F,
                             minClusterSize = minModuleSize)
table(dynamicMods)

dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

sizeGrWindow(8, 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic tree cut",
                    dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05,
                    main = "Gene dendrogram and module colors\nDistal colon")

MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

MEDiss <- 1 - cor(MEs)

METree <- hclust(as.dist(MEDiss), method = "average")

sizeGrWindow(7, 6)
pdf(paste0(plots, "distal_colon_eigengene_dendrogram.pdf"), height = 6, width = 7)
plot(METree, main = "Clustering of module eigengenes\nDistal colon",
     xlab = "", sub = "")
MEDissThres <- 0.25
abline(h = MEDissThres, col = "red")
dev.off()

merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

sizeGrWindow(12, 9)
pdf(paste0(plots, "distal_colon_cluster_dendrogram.pdf"), height = 9, width = 12)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic tree cut", "Merged dynamic"),
                    dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05)
dev.off()

moduleColors <- mergedColors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs

save(MEs, moduleLabels, moduleColors, geneTree, file = "distal_colon_WGCNA_build.RData")

message("Starting proximal colon dataset")
#PROXIMAL COLON----
datExpr <- multiExpr[[2]]$data

adjacency <- adjacency(datExpr, power = softPower, type = "signed hybrid")

TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")

sizeGrWindow(12,9)
pdf(paste0(plots, "proximal_colon_gene_clustering.pdf", height = 9, width = 12))
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity\nProximal colon",
     labels = F, hang = 0.04)
dev.off()

minModuleSize <- 30

dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = F,
                             minClusterSize = minModuleSize)
table(dynamicMods)

dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

sizeGrWindow(8, 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic tree cut",
                    dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05,
                    main = "Gene dendrogram and module colors\nProximal colon")

MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

MEDiss <- 1 - cor(MEs)

METree <- hclust(as.dist(MEDiss), method = "average")

sizeGrWindow(7, 6)
pdf(paste0(plots, "proximal_colon_eigengene_dendrogram.pdf"), height = 6, width = 7)
plot(METree, main = "Clustering of module eigengenes\nProximal colon",
     xlab = "", sub = "")
MEDissThres <- 0.25
abline(h = MEDissThres, col = "red")
dev.off()

merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

sizeGrWindow(12, 9)
pdf(paste0(plots, "proximal_colon_cluster_dendrogram.pdf"), height = 9, width = 12)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic tree cut", "Merged dynamic"),
                    dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05)
dev.off()

moduleColors <- mergedColors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs

save(MEs, moduleLabels, moduleColors, geneTree, file = "proximal_colon_WGCNA_build.RData")
