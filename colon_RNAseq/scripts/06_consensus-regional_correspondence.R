#This script performs multiple Fisher's exact tests to determine how similar consensus module assignment is to regionally-specific module assignment
#This analysis uses the consensus WGCNA from 03 and the regional WGCNAs from 04
#Cells that show a significant Bonferroni-adjusted p-value (< 0.05) mean that there is significant overlap in the gene assigment of those modules
#In other words, the regional and consensus modules that align in that cell show enough overlap to call them "similar", meaning that we can consider there to be
#good preservation of the consensus module in the individual colon segment