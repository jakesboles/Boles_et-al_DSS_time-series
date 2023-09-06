#This script conducts an analysis that evaluates how well the consensus network is preserved between both colon segments
#This analysis is wrapped neatly within the WGCNA package

library(WGCNA)
library(rstudioapi)

load(file = "colon_RNAseq/data_objects/consensus_WGCNA_input.RData")
load(file = "colon_RNAseq/data_objects/consensus_network.RData")

exprSize <- checkSets(multiExpr)
nSets <- exprSize$nSets

consMEsC <- multiSetMEs(multiExpr, universalColors = moduleColors)
MET <- consensusOrderMEs(consMEsC)

sizeGrWindow(8, 10)
par(cex = 0.9)
plotEigengeneNetworks(MET, setLabels, marDendro = c(0, 2, 2, 1),
                      marHeatmap = c(3, 3, 2, 1),
                      zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
savePlotAsImage(file = "colon_RNAseq/plots/colon_consensus_networks_preservation.png", 
                format = "png",
                height = 1000, width = 1200)