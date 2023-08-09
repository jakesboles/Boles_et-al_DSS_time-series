#This script visualizes the module preservation results in one heatmap, focusing only on modules of interest
#Interesting modules were identified in 07

library(WGCNA)
library(tidyverse)
library(ggplot2)
library(paletteer)
library(forcats)
library(pheatmap)
library(cowplot)

load("colon_RNAseq/data_objects/module_preservation.RData")

ref <- 1
test <- 2
mp <- mp.d1 
statsObs1.d <- cbind(mp$quality$observed[[ref]][[test]][, -1],
                     mp$preservation$observed[[ref]][[test]][, -1])
statsZ1.d <- cbind(mp$quality$Z[[ref]][[test]][, -1],
                   mp$preservation$Z[[ref]][[test]][, -1])
statsP1.d <- cbind(mp$preservation$log.p[[ref]][[test]][, -1])

ref <- 1
test <- 2
mp <- mp.p1 
statsObs1.p <- cbind(mp$quality$observed[[ref]][[test]][, -1],
                     mp$preservation$observed[[ref]][[test]][, -1])
statsZ1.p <- cbind(mp$quality$Z[[ref]][[test]][, -1],
                   mp$preservation$Z[[ref]][[test]][, -1])
statsP1.p <- cbind(mp$preservation$log.p[[ref]][[test]][, -1])

ref <- 1
test <- 2
mp <- mp.d2 
statsObs2.d <- cbind(mp$quality$observed[[ref]][[test]][, -1],
                     mp$preservation$observed[[ref]][[test]][, -1])
statsZ2.d <- cbind(mp$quality$Z[[ref]][[test]][, -1],
                   mp$preservation$Z[[ref]][[test]][, -1])
statsP2.d <- cbind(mp$preservation$log.p[[ref]][[test]][, -1])

ref <- 1
test <- 2
mp <- mp.p2 
statsObs2.p <- cbind(mp$quality$observed[[ref]][[test]][, -1],
                     mp$preservation$observed[[ref]][[test]][, -1])
statsZ2.p <- cbind(mp$quality$Z[[ref]][[test]][, -1],
                   mp$preservation$Z[[ref]][[test]][, -1])
statsP2.p <- cbind(mp$preservation$log.p[[ref]][[test]][, -1])

ref <- 1
test <- 2
mp <- mp.d3 
statsObs3.d <- cbind(mp$quality$observed[[ref]][[test]][, -1],
                     mp$preservation$observed[[ref]][[test]][, -1])
statsZ3.d <- cbind(mp$quality$Z[[ref]][[test]][, -1],
                   mp$preservation$Z[[ref]][[test]][, -1])
statsP3.d <- cbind(mp$preservation$log.p[[ref]][[test]][, -1])

ref <- 1
test <- 2
mp <- mp.p3 
statsObs3.p <- cbind(mp$quality$observed[[ref]][[test]][, -1],
                     mp$preservation$observed[[ref]][[test]][, -1])
statsZ3.p <- cbind(mp$quality$Z[[ref]][[test]][, -1],
                   mp$preservation$Z[[ref]][[test]][, -1])
statsP3.p <- cbind(mp$preservation$log.p[[ref]][[test]][, -1])

Zstats <- list(statsZ1.d, statsZ1.p, statsZ2.d, statsZ2.p,
               statsZ3.d, statsZ3.p)
names(Zstats) <- c("statsZ1.d", "statsZ1.p", "statsZ2.d", "statsZ2.p",
                   "statsZ3.d", "statsZ3.p")

Zsumm <- list()

for (i in seq_along(Zsumm)){
  Zsumm[[i]] <- Zstats[[i]]$Zsummary.pres
  names(Zsumm)[[i]] <- names(Zstats)[[i]]
}

colors <- rownames(statsZ1.d)

Zsumm <- list2DF(Zsumm) %>%
  mutate(module = colors) %>%
  pivot_longer(cols = !module,
               names_to = comparison,
               values_to = "stat")

#Make vectors containing module colors of interest
#These are in order in the way they appear in the figure 
immune <- c("blueviolet", "darkseagreen3", "darkorange", "orange", "orangered1",
            "red", "firebrick3")
stress <- c("coral2", "coral3", "white", "thistle2")
metab <- c("brown", "coral", "darkolivegreen4", "lightpink3", "coral1",
           "indianred4", "darkslateblue")
repair <- c("floralwhite", "lightcoral", "salmon", "tan")

#Concatenate vectors
all_modules <- c(immune, stress, metab, repair)

#Make a dataframe containing module grouping labels 
#This will be used in faceting/labeling the heatmap
moduletype <- data.frame(module = all_modules,
                         type = c(rep("Immune\nFig. 2", length(immune)),
                                  rep("Cell stress\nFig. S4", length(stress)),
                                  rep("Metabolism\nFig. S5", length(metab)),
                                  rep("Repair\nFig. S6", length(repair))))

(p <- Zsumm[Zsumm$module %in% all_modules, ] %>%
    left_join(moduletype, by = "module") %>%
    mutate(module = factor(module, 
                           levels = fct_inorder(moduletype$module))) %>%
    mutate(comparison = factor(comparison,
                               levels = c("statsZ1.d", "statsZ1.p", "statsZ2.d", 
                                          "statsZ2.p", "statsZ3.d", "statsZ3.p"),
                               labels = c("Reference:\ndistal colon\nTest:\nGSE131032",
                                          "Reference:\nproximal colon\nTest:\nGSE131032",
                                          "Reference:\ndistal colon\nTest:\nGSE168053\ndistal colon",
                                          "Reference:\nproximal colon\nTest:\nGSE168053\nproximal colon",
                                          "Reference:\ndistal colon\nTest:\nGSE210405",
                                          "Reference:\nproximal colon\nTest:\nGSE210405")),
           module = str_to_title(module)) %>%
    ggplot(aes(x = comparison, y = factor(module,
                                          levels = fct_inorder(str_to_title(moduletype$module))))) +
    geom_tile(aes(fill = stat),
              color = "white") + 
    scale_fill_stepsn(name = "Zsummary",
                      colors = c("white", "black", "#007D97", "#79D359", "#FDE333", "#D33F6A"),
                      breaks = c(2, 10, 20, 30, 40),
                      limits = c(-7, 51)) +
    scale_y_discrete(position = "right",
                     limits = rev) +
    ggtitle("Preservation of colon WGCNA consensus\nmodules in public DSS intestine datasets") +
    facet_grid(rows = vars(factor(type,
                                  levels = unique(fct_inorder(moduletype$type)))),
               scales = "free_y", 
               space = "free_y",
               switch = "y") + 
    theme_classic(base_size = 16) + 
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(color = "black", angle = 0,
                                 hjust = 0.5, vjust = 1),
      axis.text.y = element_text(color = "black"),
      axis.title = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(size = 20, hjust = 0.5)
    ))

ggsave(p,
       filename = "colon_RNAseq/plots/module_preservation_heatmap.png",
       units = "in", dpi = 600,
       height = 8, width = 11)
