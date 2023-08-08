#This script ranks genes using a fold-change shrinkage strategy via DESeq2
#Ranked gene lists will be our input for the GSPA analysis
#Fold-change shrinkage is important here and I would refer you to DESeq2's documentation for the rationale

library(biomaRt) #Bioc
library(dplyr)
library(clusterProfiler) #Bioc
library(tidyverse)
library(Orthology.eg.db) #Bioc
library(org.Mm.eg.db) #Bioc
library(org.Hs.eg.db) #Bioc
library(AnnotationDbi) #Bioc
library(DESeq2)
library(rlist)
library(apeglm)

mapfun <- function(mousegenes){
  gns <- mapIds(org.Mm.eg.db, mousegenes, "ENTREZID", "SYMBOL")
  mapped <- select(Orthology.eg.db, gns, "Homo_sapiens","Mus_musculus")
  naind <- is.na(mapped$Homo_sapiens)
  hsymb <- mapIds(org.Hs.eg.db, as.character(mapped$Homo_sapiens[!naind]), "SYMBOL", "ENTREZID")
  out <- data.frame(Mouse_symbol = mousegenes, mapped, Human_symbol = NA)
  out$Human_symbol[!naind] <- hsymb
  out
}

#INSERT DESEQ2 OBJECT LOADING HERE

c_gr1 <- as.data.frame(lfcShrink(c_dds, coef = "group_5d_DSS_vs_Untreated", type = "apeglm"))
c_gr2 <- as.data.frame(lfcShrink(c_dds, coef = "group_7d_DSS_vs_Untreated", type = "apeglm"))
c_gr3 <- as.data.frame(lfcShrink(c_dds, coef = "group_7d_DSS_2d_H2O_vs_Untreated", type = "apeglm"))
c_gr4 <- as.data.frame(lfcShrink(c_dds, coef = "group_7d_DSS_5d_H2O_vs_Untreated", type = "apeglm"))
c_gr5 <- as.data.frame(lfcShrink(c_dds, coef = "group_7d_DSS_7d_H2O_vs_Untreated", type = "apeglm"))
c_gr6 <- as.data.frame(lfcShrink(c_dds, coef = "group_7d_DSS_14d_H2O_vs_Untreated", type = "apeglm"))

m_gr1 <- as.data.frame(lfcShrink(m_dds, coef = "group_5d_DSS_vs_Untreated", type = "apeglm"))
m_gr2 <- as.data.frame(lfcShrink(m_dds, coef = "group_7d_DSS_vs_Untreated", type = "apeglm"))
m_gr3 <- as.data.frame(lfcShrink(m_dds, coef = "group_7d_DSS_2d_H2O_vs_Untreated", type = "apeglm"))
m_gr4 <- as.data.frame(lfcShrink(m_dds, coef = "group_7d_DSS_5d_H2O_vs_Untreated", type = "apeglm"))
m_gr5 <- as.data.frame(lfcShrink(m_dds, coef = "group_7d_DSS_7d_H2O_vs_Untreated", type = "apeglm"))
m_gr6 <- as.data.frame(lfcShrink(m_dds, coef = "group_7d_DSS_14d_H2O_vs_Untreated", type = "apeglm"))

h_gr1 <- as.data.frame(lfcShrink(h_dds, coef = "group_5d_DSS_vs_Untreated", type = "apeglm"))
h_gr2 <- as.data.frame(lfcShrink(h_dds, coef = "group_7d_DSS_vs_Untreated", type = "apeglm"))
h_gr3 <- as.data.frame(lfcShrink(h_dds, coef = "group_7d_DSS_2d_H2O_vs_Untreated", type = "apeglm"))
h_gr4 <- as.data.frame(lfcShrink(h_dds, coef = "group_7d_DSS_5d_H2O_vs_Untreated", type = "apeglm"))
h_gr5 <- as.data.frame(lfcShrink(h_dds, coef = "group_7d_DSS_7d_H2O_vs_Untreated", type = "apeglm"))
h_gr6 <- as.data.frame(lfcShrink(h_dds, coef = "group_7d_DSS_14d_H2O_vs_Untreated", type = "apeglm"))

s_gr1 <- as.data.frame(lfcShrink(s_dds, coef = "group_5d_DSS_vs_Untreated", type = "apeglm"))
s_gr2 <- as.data.frame(lfcShrink(s_dds, coef = "group_7d_DSS_vs_Untreated", type = "apeglm"))
s_gr3 <- as.data.frame(lfcShrink(s_dds, coef = "group_7d_DSS_2d_H2O_vs_Untreated", type = "apeglm"))
s_gr4 <- as.data.frame(lfcShrink(s_dds, coef = "group_7d_DSS_5d_H2O_vs_Untreated", type = "apeglm"))
s_gr5 <- as.data.frame(lfcShrink(s_dds, coef = "group_7d_DSS_7d_H2O_vs_Untreated", type = "apeglm"))
s_gr6 <- as.data.frame(lfcShrink(s_dds, coef = "group_7d_DSS_14d_H2O_vs_Untreated", type = "apeglm"))

gspa_genes <- read.csv("brain_RNAseq/gspa_genes.csv")
gspa_genes <- gspa_genes %>%
  pull(2)

#Create list of data frames to loop over - just reduces length of code
dfs <- list(c_gr1, c_gr2, c_gr3, c_gr4, c_gr5, c_gr6,
            s_gr1, s_gr2, s_gr3, s_gr4, s_gr5, s_gr6,
            h_gr1, h_gr2, h_gr3, h_gr4, h_gr5, h_gr6,
            m_gr1, m_gr2, m_gr3, m_gr4, m_gr5, m_gr6)
names(dfs) <- c("c_gr1", "c_gr2", "c_gr3", "c_gr4", "c_gr5", "c_gr6",
                "s_gr1", "s_gr2", "s_gr3", "s_gr4", "s_gr5", "s_gr6",
                "h_gr1", "h_gr2", "h_gr3", "h_gr4", "h_gr5", "h_gr6",
                "m_gr1", "m_gr2", "m_gr3", "m_gr4", "m_gr5", "m_gr6")
for (i in seq_along(dfs)){
  dfs[[i]] <- dfs[[i]] %>%
    rownames_to_column(var = "X") %>%
    mutate(X = convert_mouse_to_human_symbols(X)) %>% #covert mouse to human symbols
    na.omit() %>% #remove genes without human homologue 
    dplyr::select(c(X, log2FoldChange)) %>% 
    arrange(desc(log2FoldChange)) #rank by shrunken fold change
  
  dfs[[i]] <- dfs[[i]][dfs[[i]]$X %in% gspa_genes, ] #select only genes that are found in GSPA embeddings
  
  write.table(dfs[[i]], file = paste0("brain_RNAseq/gspa_input/", #write tab-separated 
                                      names(dfs)[i], ".rnk"),            #.rnmk files for export 
              row.names = F, col.names = F,
              sep = "\t",
              quote = F)
}