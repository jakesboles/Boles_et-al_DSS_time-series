#This script compares the tissue-specific network analysis performed in 08 to the consensus analysis performed in 07
#It uses Fisher's exact test (a contingency table analysis) to assess how independent gene module assignment is between networks
#We Bonferroni-adjust the p-values to reflect the number of comparisons being made
#Adjusted p-values < 0.05 mean that there is significant overlap in gene assignment, which in plain English means that
#gene assignment is similar between networks