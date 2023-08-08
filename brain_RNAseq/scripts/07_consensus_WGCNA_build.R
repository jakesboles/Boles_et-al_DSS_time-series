#This script is the meat of the analysis - building the consensus weighted co-expression analysis with the WGCNA pipeline.
#Parameters for network construction are identified in this script, the TOM is created, and genes are clustered into modules, which are then grouped as needed based on their dissimilarity.
#Note that gene clustering is done with *average linkage* 