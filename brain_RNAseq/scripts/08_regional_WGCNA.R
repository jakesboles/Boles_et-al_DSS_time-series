#This script follows the same parameters identified in 07 to create weight co-expression networks in each brain region separately. 
#The only purpose of this analysis is to later perform one way of evaluating how well the consensus network is preserved in all four regions.

library(WGCNA)

load("brain_RNAseq/consensus_WGCNA_input.RData")

plots <- "brain_RNAseq/plots/"

options(stringsAsFactors = F)
enableWGCNAThreads()

softPower <- 6

message("Starting cortex")
#Cortex WGCNA----
datExpr <- multiExpr[[1]]$data
adjacency <- adjacency(datExpr, power = softPower,
                       type = "signed hybrid")

TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")

sizeGrWindow(12, 9)
pdf(paste0(plots, "cortex_gene_clustering.pdf", height = 9, width = 12))
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity\nCortex",
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
                    main = "Gene dendrogram and module colors\nCortex")

MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

MEDiss <- 1 - cor(MEs)

METree <- hclust(as.dist(MEDiss), method = "average")

sizeGrWindow(7, 6)
pdf(paste0(plots, "cortex_eigengene_dendrogram.pdf"), height = 6, width = 7)
plot(METree, main = "Clustering of module eigengenes\nCortex",
     xlab = "", sub = "")
MEDissThres <- 0.25
abline(h = MEDissThres, col = "red")
dev.off()

merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

sizeGrWindow(12, 9)
pdf(paste0(plots, "cortex__cluster_dendrogram.pdf"), height = 9, width = 12)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic tree cut", "Merged dynamic"),
                    dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05)
dev.off()

moduleColors <- mergedColors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs

save(MEs, moduleLabels, moduleColors, geneTree, file = "cortex_WGCNA_build.RData")

message("Starting midbrain dataset")
#Midbrain WGCNA----