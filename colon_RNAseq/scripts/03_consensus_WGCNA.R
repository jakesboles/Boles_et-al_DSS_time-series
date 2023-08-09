#This script builds the consensus weighted co-expression network in our two colon segments 

library(WGCNA)

options(stringsAsFactors = F)
enableWGCNAThreads() #enable parallel processing will only work on an HPC of some sort 

plots <- "colon_RNAseq/plots/"

load("colon_RNAseq/consensus_WGCNA_input.RData")

nSets <- checkSets(multiExpr)$nSets

powers <- c(seq(4, 10, by = 1), seq(12, 20, by = 2)) 
powerTables <- vector(mode = "list", length = nSets)

for (set in 1:nSets){
  powerTables[[set]] <- list(data = pickSoftThreshold(multiExpr[[set]]$data,
                                                      powerVector = powers, 
                                                      verbose = 2,
                                                      networkType = "signed hybrid")[[2]])
}
collectGarbage()

colors <- c("palegreen4", "hotpink3")
plotCols <- c(2, 5, 6, 7)
colNames <- c("Scale Free Topology Model Fit", "Mean connectivity",
              "Median connectivity", "Max connectivity")
ylim <- matrix(NA, nrow = 2, ncol = 4)
for (set in 1:nSets){
  for (col in 1:length(plotCols)){
    ylim[1, col] <- min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = T)
    ylim[2, col] <- max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = T)
  }
}

sizeGrWindow(8, 6)
pdf(paste0(plots, "soft_power_selection.pdf"), height = 6, width = 8)
par(mfcol = c(2,2))
par(mar = c(4.2, 4.2, 2.2, 0.5))
cex1 = 0.7
for (col in 1:length(plotCols)) for (set in 1:nSets){
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20)
}
dev.off()

softPower <- 7
adjacencies <- array(0, dim = c(nSets, nGenes, nGenes))

for (set in 1:nSets) {
  adjacencies[set, , ] <- cor(multiExpr[[set]]$data, use = "p") #calculate correlation coefficient
  adjacencies[set, , ][adjacencies[set , , ] < 0] <- 0 #where coef < 0, set to 0 (hard network)
  adjacencies[set, , ] <- adjacencies[set , , ]^softPower #scale up positive coef (soft network)
} 

#Create the topological overlap matrix
TOM <- array(0, dim = c(nSets, nGenes, nGenes))
for (set in 1:nSets){
  TOM[set, , ] <- TOMsimilarity(adjacencies[set, , ], TOMType = "signed")
}

#Scaling the TOM so all 95th quantiles are equal across networks
scaleP <- 0.95
set.seed(12345)
nSamples <- as.integer(1/(1-scaleP) * 1000)
scaleSample <- sample(nGenes * (nGenes-1) / 2, size = nSamples)
TOMScalingSamples <- list()

scaleQuant <- rep(1, nSets)
scalePowers <- rep(1, nSets)

for (set in 1:nSets){
  TOMScalingSamples[[set]] <- as.dist(TOM[set, , ])[scaleSample]
  scaleQuant[set] <- quantile(TOMScalingSamples[[set]],
                              probs = scaleP, type = 8)
  if (set > 1){
    scalePowers[set] <- log(scaleQuant[1])/log(scaleQuant[set])
    TOM[set, , ] <- TOM[set, , ]^scalePowers[set]
  }
}

#plot the scaled TOMs
scaledTOMSamples <- list()
for (set in 1:nSets){
  scaledTOMSamples[[set]] <- TOMScalingSamples[[set]]^scalePowers[set]
}

sizeGrWindow(6,6)
pdf(paste0(plots, "scaled_TOMs.pdf"), height = 6, width = 6)
qqUnscaled <- qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = T,
                     cex = 0.6, xlab = paste("TOM in", setLabels[1]),
                     ylab = paste("TOM in", setLabels[2]),
                     main = "QQ plot of TOM", pch = 20)
qqScaled <- qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = F)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)
abline(a = 0, b = 1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, 
       col = c("black", "red"))
dev.off()

#Consensus is defined by the minimum weight
consensusTOM <- pmin(TOM[1, , ], TOM[2, , ], TOM[3, , ], TOM[4, , ])

#Cluster genes with average linkage
consTree <- hclust(as.dist(1-consensusTOM), method = "average")

#Gene dendrogram with a minimum cluster size of 30
minModuleSize <- 30

unmergedLabels <- cutreeDynamic(dendro = consTree, distM = 1 - consensusTOM,
                                deepSplit = 2, cutHeight = 0.995,
                                minClusterSize = minModuleSize,
                                pamRespectsDendro = F)
unmergedColors <- labels2colors(unmergedLabels)
table(unmergedLabels)

sizeGrWindow(8, 6)
plotDendroAndColors(consTree, unmergedColors, "Dynamic Tree Cut",
                    dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05)

unmergedMEs <- multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors) #calculate eigengenes
consMEDiss <- consensusMEDissimilarity(unmergedMEs) #calculate consensus dissimilarity
consMETree <- hclust(as.dist(consMEDiss), method = "average")

sizeGrWindow(7, 6)
pdf(paste0(plots, "consensus_eigengene_dendrogram.pdf"), height = 6, width = 7)
par(mfrow = c(1,1))
plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = clustering)
abline(h = 0.25, col = "red")
dev.off()

#Merge modules based on dissimilarity
#In this case, a dissimilarity of 0.25 was used as the threshold for merging
#In other words, modules that are 75% similar are moerged
merge <- mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)

moduleLabels <- merge$colors
moduleColors <- labels2colors(moduleLabels)
consMEs <- merge$newMEs

sizeGrWindow(9, 6)
pdf(paste0(plots, "consensus_cluster_dendrogram.pdf"), height = 6, width = 9)
plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors),
                    c("Unmerged", "Merged"),
                    dendroLabels = F, hang = 0.03,
                    addGuide = T, guideHang = 0.05)
dev.off()

save(consMEs, moduleColors, moduleLabels, consTree, file = "colon_RNAseq/data_objects/consensus_network.RData")

