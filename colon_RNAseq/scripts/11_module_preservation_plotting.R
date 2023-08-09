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