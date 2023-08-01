#Load version-controlled package and establish the locked environment 
#All packages were loaded and up-to-date on August 1, 2023

#Use this script to rebuild the renv in case of failures during testing/development...

library(renv)
renv::init()
renv::install(c("bioc::DESeq2", "tidyverse", "janitor", "pheatmap", "RColorBrewer",
                "ggplot2", "bioc::edgeR", "gtools", "bioc::impute", "bioc::preprocessCore",
                "bioc::GO.db", "matrixStats", "Hmisc", "splines", "foreach", 
                "doParallel", "fastcluster", "dynamicTreeCut", "survival",
                "WGCNA", "paletteer", "magrittr", "cowplot", "ggbeeswarm",
                "afex", "emmeans", "multcomp", "stringr", "scales", "ggpubr",
                "VennDiagram", "data.table", "ggh4x"))

renv::install(c("bioc::biomaRt", "bioc::org.Mm.eg.db", "bioc::org.Hs.eg.db",
                "rlist", "bioc::apeglm"))

renv::install(c("bioc::Orthology.eg.db", "bioc::clusterProfiler", "reshape2",
                "ggrepel", "bioc::vsn"))
renv::install(c("readxl", "ggpp", "psych"))

renv::install("Matrix")

renv::settings$snapshot.type("all")

renv::snapshot()
