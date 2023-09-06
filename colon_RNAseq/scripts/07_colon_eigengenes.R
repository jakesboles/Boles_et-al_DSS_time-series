#This script calculates consensus module eigengenes and plots them for all consensus modules in one large figure
#This figure was used to categorize modules according to their regulation during/after DSS 
#Ultimately, this was used to correlate brain gene signatures with colon gene signatures that we've defined temporally 
#This script will also generate a large file with gene-level information, including module membership stats and correlation information with traits from traits.csv
library(WGCNA)
library(ggplot2)
library(ggbeeswarm)
library(tidyverse)
library(rstudioapi)
library(afex)
library(emmeans)
library(multcomp)
library(stringr)

options(stringsAsFactors = F)

lnames <- load("RNAseq/r_objects3/colon_04_WGCNA_consensus_input.RData")
lnames
lnames <- load(paste0("RNAseq/r_objects3/colon_05_WGCNA_consensus_network_build_", clustering, ".RData"))
lnames

consMEsC <- multiSetMEs(multiExpr, universalColors = moduleColors)
MET <- consensusOrderMEs(consMEsC)

write.csv(consMEsC[[1]]$data, file = "colon_RNAseq/csv_outputs/distal_colon_consensus_eigengenes.csv",
          row.names = T)
write.csv(consMEsC[[2]]$data, file = "RNAseq/csv_outputs3/proximal_colon_consensus_eigengenes.csv",
          row.names = T)

me <- rbind(consMEsC[[1]]$data, consMEsC[[2]]$data)
exprSize <- checkSets(multiExpr)
exprSize
me <- me %>%
  rownames_to_column(var = "FileID")

key <- read.csv("colon_RNAseq/RNAseq_colon_metadata.csv")

df <- me %>%
  left_join(key, by = "FileID") %>%
  dplyr::select(-c(t))
names <- colnames(me)
names <- names[-1]
cols <- substring(names, 3)

df2 <- df %>%
  mutate(tissue = str_to_title(tissue) %>%
           str_replace_all("_", " ")) %>%
  pivot_longer(cols = -c(FileID,id, sex, group, tissue),
               names_to = "module", values_to = "eigengene") %>%
  mutate(module = str_remove_all(module, "ME") %>%
           str_to_title())

genes <- colnames(multiExpr[[1]]$data)

colorkey <- data.frame(module = moduleColors, gene = genes)
colorkey <- colorkey %>%
  group_by(module) %>%
  summarize(n = n())
colorkey$module <- str_to_title(colorkey$module)

df2 <- df2 %>%
  left_join(colorkey, by = "module") %>%
  tidyr::unite(module, c(module, n), sep = " (")
df2$module <- paste0(df2$module, ")")

head(df2)

ggplot(df2, aes(x = factor(group,
                           levels = c("Untreated", "5d_DSS", "7d_DSS", 
                                      "7d_DSS_2d_H2O", "7d_DSS_5d_H2O",
                                      "7d_DSS_7d_H2O", "7d_DSS_14d_H2O"),
                           labels = c("Untreated", "5d", "7d", "7d + 2d",
                                      "7d + 5d", "7d + 7d", "7d + 14d")),
                y = eigengene)) +
  geom_quasirandom(aes(fill = tissue, shape = tissue),
                   dodge.width = 1,
                   size = 3,
                   alpha = 0.8,
                   show.legend = T) +
  stat_summary(fun.y = "mean", geom = "line", linewidth = 2,
               aes(group = factor(tissue), color = tissue)) +
  facet_wrap(vars(module), ncol = 7, scales = "free_y") +
  scale_shape_manual(values = c(21, 22)) + 
  scale_fill_manual(values = c("palegreen", "lightpink")) +
  scale_color_manual(values = c("palegreen4", "hotpink3")) +
  ylab("Module eigengene") +
  theme_bw(base_size = 12) +
  theme(
    axis.title.y = element_text(size = 18), 
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(size = 16, color = "white", face = "bold"),
    legend.position = c(0.4, 0.035),
    legend.text = element_text(size = 24),
    legend.title = element_blank(),
    #legend.key.size = unit(8, "pt"),
    #axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 45, color = "black", vjust = 1, 
                               size = 14, hjust = 1),
    axis.text.y = element_text(color = "black"))
ggsave("colon_RNAseq/plots/consensus_eigengenes_over_time.png", 
       units = "in", 
       dpi = 600, 
       height = 25, width = 20)

#Eigengenes, kME, GS----
consMEs.unord <- multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey = T)
GS <- list()
kME <- list()
for (set in 1:nSets){
  GS[[set]] <- corAndPvalue(multiExpr[[set]]$data, multiTraits[[set]]$data)
  kME[[set]] <- corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data)
}

GS.metaZ <- (GS[[1]]$Z + GS[[2]]$Z)/sqrt(2)
kME.metaZ <- (kME[[1]]$Z + kME[[2]]$Z)/sqrt(2)
GS.metaP <- 2*pnorm(abs(GS.metaZ), lower.tail = F)
kME.metaP <- 2*pnorm(abs(kME.metaZ), lower.tail = F)

GSmat <- rbind(GS[[1]]$cor, GS[[2]]$cor, GS[[1]]$p, GS[[2]]$p, GS.metaZ, GS.metaP)
nTraits <- checkSets(multiTraits)$nGenes
traitNames <- colnames(multiTraits[[1]]$data)
dim(GSmat) <- c(nGenes, 6*nTraits)
rownames(GSmat) <- names(multiExpr[[1]]$data)
colnames(GSmat) <- spaste(
  c("GS.distal", "GS.proximal", "p.GS.distal", "p.GS.proximal", "Z.GS.meta", "p.GS.meta"),
  rep(traitNames, rep(6, nTraits)))

kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[1]]$p, kME[[2]]$p, kME.metaZ, kME.metaP);
MEnames = colnames(consMEs.unord[[1]]$data);
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nGenes, 6*nMEs)
rownames(kMEmat) = names(multiExpr[[1]]$data);
colnames(kMEmat) = spaste(
  c("kME.distal", "kME.proximal", "p.kME.distal", "p.kME.proximal", "Z.kME.meta", "p.kME.meta"),
  rep(MEnames, rep(6, nMEs)))

info <- data.frame(Gene = names(multiExpr[[1]]$data), ModuleLabel = moduleLabels,
                   ModuleColors = labels2colors(moduleLabels),
                   GSmat, kMEmat)
write.csv(info, file = "colon_RNAseq/csv_outputs/colon_consensus_combined_result.csv",
          row.names = F, quote = F)
