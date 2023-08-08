#This script creates many DESeq2 result objects from our DESeq2 object created in 01 and generates clean volcano plots displaying the results.

library(DESeq2)
library(tidyverse)

load("brain_RNAseq/data_objects/deseq2_objects.RData")

coldata <- colData(dds)
counts <- assay(dds)

ctx_counts <- counts[, coldata$tissue %in% "cortex"]
ctx_coldata <- coldata[coldata$FileID %in% colnames(ctx_counts), ]

hpc_counts <- counts[, coldata$tissue %in% "hippocampus"]
hpc_coldata <- coldata[coldata$FileID %in% colnames(hpc_counts), ]

mb_counts <- counts[, coldata$tissue %in% "midbrain"]
mb_coldata <- coldata[coldata$FileID %in% colnames(mb_counts), ]

str_counts <- counts[, coldata$tissue %in% "striatum"]
str_coldata <- coldata[coldata$fileID %in% colnames(str_counts), ]

ctx_dds <- DESeqDataSetFromMatrix(ctx_counts, 
                                  ctx_coldata,
                                  design = ~ group)
hpc_dds <- DESeqDataSetFromMatrix(hpc_counts,
                                  hpc_coldata,
                                  design = ~ group)
mb_dds <- DESeqDataSetFromMatrix(mb_counts,
                                 mb_coldata,
                                 design = ~ group)
str_dds <- DESeqDataSetFromMatrix(str_counts,
                                  str_coldata,
                                  design = ~ group)

ctx_dds <- DESeq(ctx_dds)
hpc_dds <- DESeq(hpc_dds)
mb_dds <- DESeq(mb_dds)
str_dds <- DESeq(str_dds)

save(ctx_dds, hpc_dds, mb_dds, str_dds,
     file = "brain_RNAseq/data_objects/regional_deseq.RData")

plots <- "brain_RNAseq/plots/"

volcano <- function(data){
  data2 <- as.data.frame(data) %>%
    mutate(expression = case_when(log2FoldChange >= 1 & padj < 0.05 ~ "Up-regulated",
                                  log2FoldChange <= -1 & padj < 0.05 ~ "Down-regulated",
                                  TRUE ~ "Unchanged"))
  
  degs <- data2 %>%
    filter(expression != "Unchanged") %>%
    summarize(n = n()) %>%
    pull(1)
  
  ggplot(data2, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(shape = 21, size = 4, aes(fill = factor(expression,
                                                       levels = c("Unchanged", "Up-regulated", "Down-regulated"))),
               alpha = 0.7) +
    scale_fill_manual(values = c("gray25", "green", "magenta")) +
    theme_bw() +
    ggtitle(paste0(degs, " DEGs")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 21),
      plot.subtitle = element_text(hjust = 0.5, size = 21),
      legend.position = "none",
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13, color = "black"),
      plot.caption = element_text(size = 16)
    )
}

data.list <- vector("list", length = 4)
for (i in seq_along(data.list)){
  data.list[[i]] <- vector("list", length = 6)
}

tissues <- c(ctx_dds, hpc_dds, mb_dds, str_dds)
groups <- c("5d_DSS", "7d_DSS", "7d_DSS_2d_H2O", "7d_DSS_5d_H2O",
            "7d_DSS_7d_H2O", "7d_DSS_14d_H2O")

for (i in seq_along(tissues)){
  for (j in seq_along(groups)){
    data.list[[i]][[j]] <- results(tissues[i],
                                   contrast = c("group", groups[j], "Untreated"))
  }
}

v <- vector("list", length = 4)
for (i in seq_along(v)){
  v[[i]] <- vector("list", length = 6)
}

for (i in 1:4){
  for (j in 1:6){
    v[[i]][[j]] <- volcano(data.list[[i]][[j]])
  }
}

ggarrange(v[[1]][[1]], v[[1]][[2]], v[[1]][[3]], v[[1]][[4]], v[[1]][[5]], v[[1]][[6]],
          v[[2]][[1]], v[[2]][[2]], v[[2]][[3]], v[[2]][[4]], v[[2]][[5]], v[[2]][[6]],
          v[[3]][[1]], v[[3]][[2]], v[[3]][[3]], v[[3]][[4]], v[[3]][[5]], v[[3]][[6]],
          v[[4]][[1]], v[[4]][[2]], v[[4]][[3]], v[[4]][[4]], v[[4]][[5]], v[[4]][[6]],
          ncol = 6, nrow = 4)
ggsave(filename = paste0(plots, "volcano_plots.png"),
       units = "in", dpi = 600,
       height = 15, width = 32)