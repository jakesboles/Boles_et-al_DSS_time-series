#This script identifies up- and down-regulated genes from each brain region at each time point and creates simple Venn diagrams for them
#Genes must have adjusted p-value < 0.05 and |log2FoldChange| > 1 to be included here. 

library(tidyverse)
library(VennDiagram)
library(DESeq2)
library(rlist)

load("brain_RNAseq/data_objects/regional_deseq.RData")

data.list <- vector("list", length = 4)
for (i in seq_along(data.list)){
  data.list[[i]] <- vector("list", length = 6)
}

dfs <- list(ctx_dds, hpc_dds, mb_dds, str_dds)
groups <- c("5d_DSS", "7d_DSS", "7d_DSS_2d_H2O", "7d_DSS_5d_H2O",
            "7d_DSS_7d_H2O", "7d_DSS_14d_H2O")

tissue <- c("Cortex", "Hippocampus", "Midbrain", "Striatum")
group_num <- c(1, 2, 3, 4, 5, 6)

for (i in seq_along(dfs)){
  for (j in seq_along(groups)){
    data.list[[i]][[j]] <- results(dfs[[i]],
                                   contrast = c("group", groups[j], "Untreated"))
    
    data.list[[i]][[j]] <- data.list[[i]][[j]] %>%
      as.data.frame() %>%
      rownames_to_column(var = "X") %>%
      mutate(group = group_num[j] %>%
               as.numeric(),
             tissue = tissue[i]) %>%
      filter(abs(log2FoldChange) > 1 & padj < 0.05)
  }
}

vd <- function(direction = c("up", "down"), .group){
  direction_title <- str_to_title(direction)
  
  if(direction == "up") {
    ctx <- data.list[[1]][[.group]] %>% filter(log2FoldChange > 0) %>% pull(X)
    hpc <- data.list[[2]][[.group]] %>% filter(log2FoldChange > 0) %>% pull(X)
    mb <- data.list[[3]][[.group]] %>% filter(log2FoldChange > 0) %>% pull(X)
    str <- data.list[[4]][[.group]] %>% filter(log2FoldChange > 0) %>% pull(X)
  }
  if(direction == "down") {
    ctx <- data.list[[1]][[.group]] %>% filter(log2FoldChange < 0) %>% pull(X)
    hpc <- data.list[[2]][[.group]] %>% filter(log2FoldChange < 0) %>% pull(X)
    mb <- data.list[[3]][[.group]] %>% filter(log2FoldChange < 0) %>% pull(X)
    str <- data.list[[4]][[.group]] %>% filter(log2FoldChange < 0) %>% pull(X)
  }
  
  venn.diagram(
    x = list(ctx, hpc, mb, str),
    category.names = c("Cortex", "Hippo-\ncampus", "Midbrain", "Striatum"),
    filename = paste0("brain_RNAseq/plots/venn_degs_group_", .group, "_", direction, ".png"),
    output = T,
    imagetype = "png",
    units = "in",
    height = 2, width = 4,
    resolution = 600,
    lwd = 2,
    lty = 'blank',
    col = c("#4E79A7", "#CCB22B", "#59A14F", "#E15759"),
    fill = c("#A0CBE8", "#F0E442", "#8CD17D", "#FF9D9A"),
    fontfamily = "Helvetica",
    cex = 1.5,
    fontface = "bold",
    cat.cex = 1,
    cat.fontfamily = "Helvetica",
    disable.logging = T,
    main = paste0(direction_title, "-regulated genes"),
    main.fontfamily = "Helvetica",
    main.cex = 1.5,
    print.mode = c("raw")
  )
  
  max_length <- max(c(length(ctx), length(hpc), length(mb), length(str)))
  
  list <- data.frame(cortex = c(ctx, rep(NA, max_length - length(ctx))),
                     hippocampus = c(hpc, rep(NA, max_length - length(hpc))),
                     midbrain = c(mb, rep(NA, max_length - length(mb))),
                     striatum = c(str, rep(NA, max_length - length(str))))
  
  print(list)
  write.csv(list, file = paste0("brain_RNAseq/csv_outputs/", direction, "_degs_group_", .group, ".csv"),
            row.names = T)
}

for (i in 1:6){
  vd("up", i)
  vd("down", i)
}
