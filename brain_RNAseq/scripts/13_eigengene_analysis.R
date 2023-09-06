#This script will calculate and plot consensus module eigengenes by group/tissue
#Here, individual eigengenes are plotted from the modules of interest we've highlighted in 11
#We also generate a full faceted figure with every consensus module the WGCNA revealed
#and generate a large .csv file with gene-level information including kME, gene-trait significance, and p-values
library(WGCNA)
library(ggplot2)
library(ggbeeswarm)
library(tidyverse)
library(rstudioapi)
library(afex)
library(emmeans)
library(multcomp)
library(stringr)
library(scales)
library(ggpubr)

load(file = "brain_RNAseq/data_objects/consensus_WGCNA_input.RData")
load(file = "brain_RNAseq/data_objects/consensus_network.RData")

options(stringsAsFactors = F)

consMEsc <- multiSetMEs(multiExpr, universalColors = moduleColors)
MET <- consensusOrderMEs(consMEsC)

write.csv(consMEsC[[1]]$data, file = "RNAseq/csv_outputs3/cortex_consensus_eigengenes.csv",
          row.names = T)
write.csv(consMEsC[[2]]$data, file = "RNAseq/csv_outputs3/midbrain_consensus_eigengenes.csv",
          row.names = T)
write.csv(consMEsC[[3]]$data, file = "RNAseq/csv_outputs3/hippocampus_consensus_eigengenes.csv",
          row.names = T)
write.csv(consMEsC[[4]]$data, file = "RNAseq/csv_outputs3/striatum_consensus_eigengenes.csv",
          row.names = T)

me <- rbind(consMEsC[[1]]$data, consMEsC[[2]]$data, consMEsC[[3]]$data, consMEsC[[4]]$data)
me <- me %>%
  rownames_to_column(var = "FileID")

key <- read.csv("brain_RNAseq/RNAseq_brain_metadata.csv")

df <- me %>%
  left_join(key, by = "FileID") %>%
  dplyr::select(-c(t))
names <- colnames(me)
names <- names[-1]
cols <- substring(names, 3)

df2 <- df %>%
  mutate(tissue = str_to_title(tissue)) %>%
  pivot_longer(cols = -c(FileID,id, sex, group, tissue),
               names_to = "module", values_to = "eigengene") %>%
  mutate(module = str_remove_all(module, "ME") %>%
           str_to_title())

#Plotting all eigengenes in a big facetted ggplot:
ggplot(df2, aes(x = factor(group,
                           levels = c("Untreated", "5d_DSS", "7d_DSS", 
                                      "7d_DSS_2d_H2O", "7d_DSS_5d_H2O",
                                      "7d_DSS_7d_H2O", "7d_DSS_14d_H2O"),
                           labels = c("Untreated", "5d", "7d", "7d + 2d",
                                      "7d + 5d", "7d + 7d", "7d + 14d")), 
                y = eigengene)) +
  geom_quasirandom(aes(fill = tissue),
                   dodge.width = 1, shape = 21,
                   show.legend = F) +
  stat_summary(fun.y = "mean", geom = "line", linewidth = 2,
               aes(group = factor(tissue), color = tissue)) +
  facet_wrap(vars(module), nrow = 4, scales = "free_y") +
  scale_fill_manual(values = c("#A0CBE8", "#F0E442", "#8CD17D", "#FF9D9A")) +
  scale_color_manual(values = c("#4E79A7", "#CCB22B", "#59A14F", "#E15759")) +
  ylab("Eigengene") +
  theme_bw(base_size = 12) +
  theme(
    axis.title.y = element_text(size = 18), 
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(size = 12, color = "white", face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 18),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank())
ggsave("brain_RNAseq/plots/consensus_eigengenes_over_time.png", 
       units = "in", 
       dpi = 600, 
       height = 10, width = 15)

#Write a function to plot particular eigengenes of interest based on findings from 11 and 12
me_plot <- function(color, titlecolor = NULL){
  Module <- paste0("ME", color)
  if(is.null(titlecolor)) titlecolor <- color
  
  df2 %>%
    subset(module == Module) %>%
    ggplot(aes(x = factor(group,
                          levels = c("Untreated", "5d_DSS", "7d_DSS", "7d_DSS_2d_H2O",
                                     "7d_DSS_5d_H2O", "7d_DSS_7d_H2O", "7d_DSS_14d_H2O"),
                          labels = c("Untreated", "5d", "7d", "7d + 2d",
                                     "7d + 5d", "7d + 7d", "7d + 14d")), 
               y = eigengene)) +
    geom_quasirandom(aes(fill = tissue, shape = tissue),
                     dodge.width = 0.5,
                     show.legend = T, size = 5) +
    stat_summary(fun.y = "mean", geom = "line", linewidth = 2,
                 aes(group = factor(tissue), color = tissue),
                 position = position_dodge(width = 0.5)) +
    scale_fill_manual(values = c("#A0CBE8", "#F0E442", "#8CD17D", "#FF9D9A")) +
    scale_color_manual(values = c("#4E79A7", "#CCB22B", "#59A14F", "#E15759")) +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    scale_x_discrete(labels = label_wrap(9)) +
    ylab("Eigengene") +
    ggtitle(paste0(str_to_title(color), " module expression")) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 21, color = titlecolor, face = "bold"),
      axis.title.y = element_text(size = 18), 
      strip.background = element_rect(fill = "black"),
      strip.text = element_text(size = 12, color = "white", face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 16),
      legend.title = element_blank(),
      axis.text.x = element_text(size = 14, color = "black",
                                 angle = 24, hjust = 1, vjust = 1),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 12)
    )
}

me_cyan <- me_plot("cyan", "cyan3")
me_paleturquoise <- me_plot("paleturquoise", "paleturquoise4")

me_darkgrey <- me_plot("darkgrey")

me_up <- ggarrange(me_cyan, me_darkgrey, me_paleturquoise,
                   common.legend = T,
                   legend = "right",
                   ncol = 3, nrow = 1)
me_up
ggsave("brain_RNAseq/plots/consensus_up_me.png",
       units = "in", dpi = 600,
       height = 6, width = 18)

me_plot("midnightblue")
ggsave("brain_RNAseq/plots/consensus_midnightblue_me.png",
       units = "in", dpi = 600,
       height = 6, width = 6)

#kME, GS, and p-values----
consMEs.unord <- multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey = T)
GS <- list()
kME <- list()
for (set in 1:nSets){
  GS[[set]] <- corAndPvalue(multiExpr[[set]]$data, multiTraits[[set]]$data)
  kME[[set]] <- corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data)
}

GS.metaZ <- (GS[[1]]$Z + GS[[2]]$Z + GS[[3]]$Z + GS[[4]]$Z)/sqrt(4)
kME.metaZ <- (kME[[1]]$Z + kME[[2]]$Z + kME[[3]]$Z + kME[[4]]$Z)/sqrt(4)
GS.metaP <- 2*pnorm(abs(GS.metaZ), lower.tail = F)
kME.metaP <- 2*pnorm(abs(kME.metaZ), lower.tail = F)

GSmat <- rbind(GS[[1]]$cor, GS[[2]]$cor, GS[[3]]$cor, GS[[4]]$cor, GS[[1]]$p, GS[[2]]$p, GS[[3]]$p, GS[[4]]$p, GS.metaZ, GS.metaP)
nTraits <- checkSets(multiTraits)$nGenes
traitNames <- colnames(multiTraits[[1]]$data)
dim(GSmat) <- c(nGenes, 10*nTraits)
rownames(GSmat) <- names(multiExpr[[1]]$data)
colnames(GSmat) <- spaste(
  c("GS.cortex", "GS.midbrain", "GS.hippocampus", "GS.striatum", "p.GS.cortex", "p.GS.midbrain", "p.GS.hippocampus", "p.GS.striatum", "Z.GS.meta", "p.GS.meta"),
  rep(traitNames, rep(10, nTraits)))

kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[3]]$cor, kME[[4]]$cor, kME[[1]]$p, kME[[2]]$p, kME[[3]]$p, kME[[4]]$p, kME.metaZ, kME.metaP);
MEnames = colnames(consMEs.unord[[1]]$data);
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nGenes, 10*nMEs)
rownames(kMEmat) = names(multiExpr[[1]]$data);
colnames(kMEmat) = spaste(
  c("kME.cortex", "kME.midbrain", "kME.hippocampus", "kME.striatum", "p.kME.cortex", "p.kME.midbrain", "p.kME.hippocampus", "p.kME.striatum", "Z.kME.meta", "p.kME.meta"),
  rep(MEnames, rep(10, nMEs)))

info <- data.frame(Gene = names(multiExpr[[1]]$data), ModuleLabel = moduleLabels,
                   ModuleColors = labels2colors(moduleLabels),
                   GSmat, kMEmat)
write.csv(info, file = paste0("RNAseq/csv_outputs3/brain_consensus_combined_results_", clustering, ".csv"),
          row.names = F, quote = F)
