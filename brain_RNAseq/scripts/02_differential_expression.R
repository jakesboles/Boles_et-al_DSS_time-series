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