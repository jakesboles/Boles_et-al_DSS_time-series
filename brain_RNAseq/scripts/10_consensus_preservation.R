#This script calculates module preservation statistics for the consensus modules found in the WGCNA
#Effectively, we're evaluating how similar the modules are between one network (i.e., tissue) and another
#This analysis is wrapped neatly in the WGCNA package
library(WGCNA)
library(rstudioapi)

load(file = "brain_RNAseq/data_objects/consensus_WGCNA_input.RData")
load(file = "brain_RNAseq/data_objects/consensus_network.RData")

exprSize <- checkSets(multiExpr)
nSets <- exprSize$nSets

consMEsC <- multiSetMEs(multiExpr, universalColors = moduleColors)
MET <- consensusOrderMEs(consMEsC)

sizeGrWindow(8, 10)
par(cex = 0.9)
plotEigengeneNetworks(MET, setLabels, marDendro = c(0, 2, 2, 1),
                      marHeatmap = c(3, 3, 2, 1),
                      zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
savePlotAsImage(file = "brain_RNAseq/plots/brain_consensus_networks_preservation.png", 
                format = "png",
                height = 1000, width = 1200)