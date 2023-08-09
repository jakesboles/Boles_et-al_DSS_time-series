#This script visualizes the module preservation results in one heatmap, focusing only on modules of interest
#Interesting modules were identified in 07

library(WGCNA)
library(tidyverse)
library(ggplot2)
library(paletteer)
library(forcats)
library(pheatmap)
library(ggh4x)
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
