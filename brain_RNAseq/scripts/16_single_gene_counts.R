#This will use the DESeq2 object created in 01 to plot the normalized counts of individual genes
#In the paper, genes of interest were identified based on module membership and functional pathway attribution according to anRichment (13)

library(DESeq2)
library(WGCNA)
library(ggplot2)
library(ggbeeswarm)
library(rlist)
library(afex)
library(multcomp)
library(emmeans)
library(tidyverse)
library(scales)
library(ggpubr)

load("brain_RNAseq/data_objects/deseq2_objects.RData")
traits <- read.csv("traits.csv")

temp <- data.frame(id = rep(unique(traits$id), 4),
                   tissue = c(rep("Cortex", 69),
                              rep("Hippocampus", 69),
                              rep("Midbrain", 69),
                              rep("Striatum", 69))) %>%
  left_join(traits, by = "id") %>%
  dplyr::select(c(id, group, tissue)) %>%
  mutate(group = factor(group,
                        levels = c("Untreated", "5d_DSS", "7d_DSS", "7d_DSS_2d_H2O",
                                   "7d_DSS_5d_H2O", "7d_DSS_7d_H2O", "7d_DSS_14d_H2O")))

plot_counts <- function(gene){ #this one will make the huge panel in Supp Fig 13
  d <- plotCounts(dds, gene = gene, intgroup = c("tissue", "group"),
                  returnData = T)
  d$tissue <- factor(str_to_title(d$tissue))
  
  d <- d %>%
    rownames_to_column(var = "FileID") %>%
    left_join(brain_key, by = "FileID") %>%
    right_join(temp, by = c("tissue", "id", "group")) %>%
    mutate(id = factor(id))
  
  tissue.only.lm <- lmer(count ~ tissue + (1 | id),
                         data = d,
                         na.action = na.exclude,
                         REML = F)
  maineff.lm <- lmer(count ~ group + tissue + (1 | id),
                     data = d,
                     na.action = na.exclude,
                     REML = F)
  full.lm <- lmer(count ~ group * tissue + (1 | id),
                  data = d,
                  na.action = na.exclude,
                  REML = F)
  message("Significant effect of group?")
  print(anova(tissue.only.lm, maineff.lm))
  message("Significant interaction between tissue and group?")
  print(anova(maineff.lm, full.lm))
  
  as.data.frame(anova(tissue.only.lm, maineff.lm))[2,] %>%
    pull(8) %>%
    round(digits = 3) -> pval1
  as.data.frame(anova(tissue.only.lm, maineff.lm))[2,] %>%
    pull(7) -> df1
  as.data.frame(anova(tissue.only.lm, maineff.lm))[2,] %>%
    pull(6) %>%
    round(digits = 3) -> chisq1
  
  as.data.frame(anova(maineff.lm, full.lm))[2, ] %>%
    pull(8) %>%
    round(digits = 3) -> pval2
  as.data.frame(anova(maineff.lm, full.lm))[2, ] %>%
    pull(7) -> df2
  as.data.frame(anova(maineff.lm, full.lm))[2, ] %>%
    pull(6) %>%
    round(digits = 3) -> chisq2
  
  pval1 <- if_else(pval1 < 0.001, "< 0.001", paste0("= ", pval1))
  pval2 <- if_else(pval2 < 0.001, " < 0.001", paste0("= ", pval2))
  
  emm <- emmeans(full.lm, specs = pairwise ~ group | tissue)
  cld <- cld(emm, Letters = letters) %>%
    mutate(.group = str_remove_all(.group, " ")) %>%
    dplyr::rename("letters" = ".group")
  
  y_add <- d %>%
    group_by(group, tissue) %>%
    summarize(x = sd(count, na.rm = T)) %>%
    arrange(desc(x)) %>%
    .[3, 3] %>%
    as.numeric()
  
  cld.y <- d %>%
    group_by(group, tissue) %>%
    summarize(y = max(count, na.rm = T) + y_add) %>%
    left_join(cld, by = c("group", "tissue"))
  
  d %>%
    ggplot(aes(y = count,
               x = factor(group,
                          levels = c("Untreated", "5d_DSS", "7d_DSS", "7d_DSS_2d_H2O",
                                     "7d_DSS_5d_H2O", "7d_DSS_7d_H2O", "7d_DSS_14d_H2O"),
                          labels = c("Untreated", "5d", "7d", "7d + 2d",
                                     "7d + 5d", "7d + 7d", "7d + 14d")))) +
    geom_quasirandom(aes(color = tissue, shape = tissue, fill = tissue),
                     dodge.width = 0.8,
                     size = 4) +
    stat_summary(aes(color = tissue),
                 fun = "mean", geom = "crossbar",
                 width = 0.5, 
                 position = position_dodge(width = 0.8),
                 show.legend = F) +
    stat_summary(aes(color = tissue),
                 fun.data = "mean_se", geom = "errorbar",
                 width = 0.5, 
                 linewidth = 1.2,
                 position = position_dodge(width = 0.8),
                 show.legend = F) +
    scale_fill_manual(values = c("#A0CBE8", "#F0E442", "#8CD17D", "#FF9D9A")) +
    scale_color_manual(values = c("#4E79A7", "#CCB22B", "#59A14F", "#E15759")) +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    facet_wrap(vars(tissue), nrow = 1) + 
    scale_x_discrete(labels = label_wrap(9)) +
    ggtitle(gene) +
    labs(y = "Normalized counts +\n1/2 pseudocount",
         caption = 
           bquote(atop("Effect of group:" ~ chi * "[" * .(df1) * "] =" ~ .(chisq1) * ";" ~ 
                         italic("p") ~ .(pval1), "Tissue x group interaction:" ~ 
                         chi * "[" * .(df2) * "] =" ~ .(chisq2) * ";" ~ 
                         italic("p") ~ .(pval2)))) +
    theme_bw(base_size = 12) +
    theme(
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      plot.title = element_text(size = 20, face = "bold.italic", hjust = 0.5),
      axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 16),
      legend.text = element_text(size = 16),
      strip.background = element_rect(fill = "black"),
      strip.text = element_text(color = "white", size = 14, face = "bold"),
      legend.position = "none",
      plot.caption = element_text(hjust = 0.5, size = 16, face = "bold")
    ) -> p
  
  p <- p +
    geom_text(data = cld.y, aes(label = letters,
                                y = y,
                                color = tissue),
              position = position_dodge(width = 0.8),
              show.legend = F,
              size = 6, face = "bold")
  
  return(p) 
}

ggarrange(plot_counts("Ada"), plot_counts("Alas2"), plot_counts("Bhlhe40"),
          plot_counts("C4a"), plot_counts("Cd14"), plot_counts("Ch25h"),
          plot_counts("Ctla2a"), plot_counts("Dusp5"), plot_counts("Ezr"),
          plot_counts("Fbln5"), plot_counts("Fcgr4"), plot_counts("Fosb"),
          plot_counts("Gimap5"), plot_counts("Hba-a1"), plot_counts("Hbb-bt"), 
          plot_counts("Icosl"), plot_counts("Il4ra"), plot_counts("Isg20"), 
          plot_counts("Itgad"), plot_counts("Jak3"), plot_counts("Lrg1"), 
          plot_counts("Lyve1"), plot_counts("Maff"), plot_counts("Mcl1"), 
          plot_counts("Mmp8"), plot_counts("Nrros"), plot_counts("Pecam1"), 
          plot_counts("S100a9"), plot_counts("Sat1"), plot_counts("Sdc4"), 
          plot_counts("Selp"), plot_counts("Sgk1"), plot_counts("Slc4a1"), 
          plot_counts("Spsb1"), plot_counts("Trim10"),
          ncol = 4, nrow = 9)
ggsave("brain_RNAseq/plots/counts_immune.png",
       units = "in", dpi = 600,
       height = 35, width = 40)

plot_counts <- function(gene){ #this one switches the formatting slightly to get the figures in Fig 6 and 7
  d <- plotCounts(dds, gene = gene, intgroup = c("tissue", "group"),
                  returnData = T)
  d$tissue <- factor(str_to_title(d$tissue))
  
  d <- d %>%
    rownames_to_column(var = "FileID") %>%
    left_join(brain_key, by = "FileID") %>%
    right_join(temp, by = c("tissue", "id", "group")) %>%
    mutate(id = factor(id))
  
  tissue.only.lm <- lmer(count ~ tissue + (1 | id),
                         data = d,
                         na.action = na.exclude,
                         REML = F)
  maineff.lm <- lmer(count ~ group + tissue + (1 | id),
                     data = d,
                     na.action = na.exclude,
                     REML = F)
  full.lm <- lmer(count ~ group * tissue + (1 | id),
                  data = d,
                  na.action = na.exclude,
                  REML = F)
  message("Significant effect of group?")
  print(anova(tissue.only.lm, maineff.lm))
  message("Significant interaction between tissue and group?")
  print(anova(maineff.lm, full.lm))
  
  as.data.frame(anova(tissue.only.lm, maineff.lm))[2,] %>%
    pull(8) %>%
    round(digits = 3) -> pval1
  as.data.frame(anova(tissue.only.lm, maineff.lm))[2,] %>%
    pull(7) -> df1
  as.data.frame(anova(tissue.only.lm, maineff.lm))[2,] %>%
    pull(6) %>%
    round(digits = 3) -> chisq1
  
  as.data.frame(anova(maineff.lm, full.lm))[2, ] %>%
    pull(8) %>%
    round(digits = 3) -> pval2
  as.data.frame(anova(maineff.lm, full.lm))[2, ] %>%
    pull(7) -> df2
  as.data.frame(anova(maineff.lm, full.lm))[2, ] %>%
    pull(6) %>%
    round(digits = 3) -> chisq2
  
  pval1 <- if_else(pval1 < 0.001, "< 0.001", paste0("= ", pval1))
  pval2 <- if_else(pval2 < 0.001, " < 0.001", paste0("= ", pval2))
  
  emm <- emmeans(full.lm, specs = pairwise ~ group | tissue)
  cld <- cld(emm, Letters = letters) %>%
    mutate(.group = str_remove_all(.group, " ")) %>%
    dplyr::rename("letters" = ".group")
  
  y_add <- d %>%
    group_by(group, tissue) %>%
    summarize(x = sd(count, na.rm = T)) %>%
    arrange(desc(x)) %>%
    .[3, 3] %>%
    as.numeric()
  
  cld.y <- d %>%
    group_by(group, tissue) %>%
    summarize(y = max(count, na.rm = T) + y_add) %>%
    left_join(cld, by = c("group", "tissue"))
  
  d %>%
    ggplot(aes(y = count,
               x = factor(group,
                          levels = c("Untreated", "5d_DSS", "7d_DSS", "7d_DSS_2d_H2O",
                                     "7d_DSS_5d_H2O", "7d_DSS_7d_H2O", "7d_DSS_14d_H2O"),
                          labels = c("Untreated", "5d", "7d", "7d + 2d",
                                     "7d + 5d", "7d + 7d", "7d + 14d")))) +
    geom_quasirandom(aes(color = tissue, shape = tissue, fill = tissue),
                     dodge.width = 0.8,
                     size = 4) +
    stat_summary(aes(color = tissue),
                 fun = "mean", geom = "crossbar",
                 width = 0.5, 
                 position = position_dodge(width = 0.8),
                 show.legend = F) +
    stat_summary(aes(color = tissue),
                 fun.data = "mean_se", geom = "errorbar",
                 width = 0.5, 
                 linewidth = 1.2,
                 position = position_dodge(width = 0.8),
                 show.legend = F) +
    scale_fill_manual(values = c("#A0CBE8", "#F0E442", "#8CD17D", "#FF9D9A")) +
    scale_color_manual(values = c("#4E79A7", "#CCB22B", "#59A14F", "#E15759")) +
    scale_shape_manual(values = c(21, 22, 23, 24)) +
    facet_wrap(vars(tissue), nrow = 2) + 
    scale_x_discrete(labels = label_wrap(9)) +
    ggtitle(gene) +
    labs(y = "Normalized counts +\n1/2 pseudocount",
         caption = 
           bquote(atop("Effect of group:" ~ chi * "[" * .(df1) * "] =" ~ .(chisq1) * ";" ~ 
                         italic("p") ~ .(pval1), "Tissue x group:" ~ 
                         chi * "[" * .(df2) * "] =" ~ .(chisq2) * ";" ~ 
                         italic("p") ~ .(pval2)))) +
    theme_bw(base_size = 12) +
    theme(
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      plot.title = element_text(size = 20, face = "bold.italic", hjust = 0.5),
      axis.text.x = element_text(size = 14, color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 16),
      legend.text = element_text(size = 16),
      strip.background = element_rect(fill = "black"),
      strip.text = element_text(color = "white", size = 14, face = "bold"),
      legend.position = "none",
      plot.caption = element_text(hjust = 0.5, size = 14, face = "bold")
    ) -> p
  
  p <- p +
    geom_text(data = cld.y, aes(label = letters,
                                y = y,
                                color = tissue),
              show.legend = F,
              size = 5)
  
  return(p) 
}

fig6_save <- function(gene){ #function to quickly save plots for Fig 6 with appropriate dimensions
  ggsave(paste0("brain_RNAseq/plots/counts_", gene,
                ".png"),
         units = "in", dpi = 600,
         height = 7, width = 5)
}
plot_counts("Nfkbia")
fig6_save("nfkbia")

plot_counts("Fas")
fig6_save("fas")

plot_counts("Il12rb1")
fig6_save("il12rb1")

plot_counts("Nfe2l2")
fig6_save("nfe2l2")

plot_counts("Hbb-bs")
fig6_save("hbbbs")

fig7_save <- function(gene){ #function to quickly save plots for Fig 7 with appropriate dimensions
  ggsave(paste0("brain_RNAseq/plots/counts_", gene,
                ".png"),
         units = "in", dpi = 600,
         height = 6, width = 6)
}

(mog <- plot_counts("Mog"))
fig7_save("mog")

(aspa <- plot_counts("Aspa"))
fig7_save("aspa")

(olig1 <- plot_counts("Olig1"))
fig7_save("Olig1")

(mbp <- plot_counts("Mbp")) 
fig7_save("mbp")

(cnp <- plot_counts("Cnp"))
fig7_save("cnp")

(mag <- plot_counts("Mag"))
fig7_save("mag")

(sox10 <- plot_counts("Sox10"))
fig7_save("sox10")

(opalin <- plot_counts("Opalin")) 
fig7_save("opalin")

(gjb1 <- plot_counts("Gjb1"))
fig7_save("gjb1")

(gjc2 <- plot_counts("Gjc2"))
fig7_save("gjc2")

