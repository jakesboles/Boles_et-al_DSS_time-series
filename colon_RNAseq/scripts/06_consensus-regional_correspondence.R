#This script performs multiple Fisher's exact tests to determine how similar consensus module assignment is to regionally-specific module assignment
#This analysis uses the consensus WGCNA from 03 and the regional WGCNAs from 04
#Cells that show a significant Bonferroni-adjusted p-value (< 0.05) mean that there is significant overlap in the gene assigment of those modules
#In other words, the regional and consensus modules that align in that cell show enough overlap to call them "similar", meaning that we can consider there to be
#good preservation of the consensus module in the individual colon segment

library(WGCNA)
library(tidyverse)
library(ggplot2)
library(reshape2)

options(stringsAsFactors = F)

tissues <- c("distal_colon", "proximal_colon")

for (i in seq_along(tissues)){
  load(file = paste0("colon_RNAseq/data_objects/", tissues[i], "_WGCNA_build.RData"))
  
  spLabels <- moduleLabels
  spColors <- moduleColors
  spTree <- geneTree
  spMEs <- orderMEs(MEs, greyName = "MEO")
  
  load(file = "colon_RNAseq/data_objects/consensus_WGCNA_input.RData")
  load(file = "colon_RNAseq/data_objects/consensus_network.RData")
  
  spModuleLabels <- substring(names(spMEs), 3)
  consModuleLabels <- substring(names(consMEs[[1]]$data), 3)
  
  spModules <- spModuleLabels
  consModules <- labels2colors(as.numeric(consModuleLabels))
  
  NspMods <- length(spModules)
  NconsMods <- length(consModules)
  
  pTable <- matrix(0, nrow = NspMods, ncol = NconsMods)
  countTbl <- matrix(0, nrow = NspMods, ncol = NconsMods)
  
  bonf <- NspMods * NconsMods
  
  for (dmod in 1:NspMods){
    for (cmod in 1:NconsMods){
      spMembers <- (spColors == spModules[dmod])
      consMembers <- (moduleColors == consModules[cmod])
      pTable[dmod, cmod] = bonf * (fisher.test(spMembers, consMembers, alternative = "greater")$p.value)
      countTbl[dmod, cmod] = sum(spColors == spModules[dmod] & moduleColors == consModules[cmod])
    }
  }
  
  spModTotals <- apply(countTbl, 1, sum)
  consModTotals <- apply(countTbl, 2, sum)
  
  pTable[pTable < 0.05] <- "P < 0.05"
  pTable[pTable != 0.05] <- "P > 0.05"
  
  consL <- length(consModules)
  spL <- length(spModules)
  
  pTable_melt <- pTable %>%
    as.data.frame()
  rownames(pTable_melt) <- spModuleLabels
  colnames(pTable_melt) <- consModules
  
  countTbl_melt <- countTbl %>%
    as.data.frame() %>%
    rownames_to_column(var = "y") %>%
    mutate(y = as.numeric(y)) %>%
    mutate(y = (spL + 1) - y) %>%
    pivot_longer(2:(ncol(countTbl) + 1),
                 names_to = "x")
  countTbl_melt$x <- str_replace_all(countTbl_melt$x, "V", "")
  countTbl_melt$x <- as.numeric(countTbl_melt$x)
  
  pTable_melt <- pTable_melt %>%
    rownames_to_column(var = "spMod") %>%
    pivot_longer(2:(ncol(pTable_melt) + 1),
                 names_to = "consMod")
  pTable_melt$spMod <- factor(pTable_melt$spMod,
                              levels = rev(unique(pTable_melt$spMod)))
  pTable_melt$consMod <- factor(pTable_melt$consMod,
                                levels = unique(pTable_melt$consMod))
  
  spGenes <- data.frame(spModules, spModTotals)
  consGenes <- data.frame(consModules, consModTotals)
  
  dev.off()
  
  pTable_melt %>%
    dplyr::rename("consModules" = "consMod", "spModules" = "spMod") %>%
    left_join(spGenes, by = "spModules") %>%
    left_join(consGenes, by = "consModules") %>%
    mutate(spModTotals = as.character(spModTotals)) %>%
    mutate(consModTotals = as.character(consModTotals)) %>%
    unite(col = "spModules", c(mbModules, spModTotals), sep = ": ") %>%
    unite(col = "consModules", c(consModules, consModTotals), sep = ": ") %>%
    mutate(spModules = factor(spModules,
                              levels = rev(unique(spModules)))) %>%
    mutate(consModules = factor(consModules,
                                levels = unique(consModules))) %>%
    ggplot(aes(consModules,spbModules)) + 
    geom_tile(aes(fill = value), color = "black") +
    scale_fill_manual(values = c("red", "white")) + 
    labs(y = paste0(str_to_title(tissues[i]), "-specific modules"), x = "Consensus modules") +
    ggtitle(paste0("Correspondence of ", tissues[i], "-specific modules and colon consensus modules")) +
    geom_text(data = countTbl_melt, aes(x = x, y = y, label = value),
              size = 3) +
    theme_classic() +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(color = "black")
    )
  
  ggsave(paste0("colon_RNAseq/plots/consensus-", tissues[i], "_correspondence.png"),
         units = "in", dpi = 600,
         height = 8, width = 10)
}