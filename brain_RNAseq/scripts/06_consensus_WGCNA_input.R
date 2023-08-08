#This script builds the input for the consensus WGCNA, and creates a layered trait object for later use. 

library(WGCNA)
library(tidyverse)
library(DESeq2)

load("brain_RNAseq/data_objects/deseq2_objects.RData")

#Creating the layered expression data object for input
#This should look familiar from 01
ctx <- vsd[, colnames(vsd) %in% brain_key[grep("cortex", brain_key$tissue), ]$FileID]
mb <- vsd[, colnames(vsd) %in% brain_key[grep("midbrain", brain_key$tissue), ]$FileID]
hpc <- vsd[, colnames(vsd) %in% brain_key[grep("hippocampus", brain_key$tissue), ]$FileID]
str <- vsd[, colnames(vsd) %in% brain_key[grep("striatum", brain_key$tissue), ]$FileID]  

nSets <- 4
setLabels <- c("Cortex", "Midbrain", "Hippocampus", "Striatum")
shortLabels <- c("CTX", "MB", "HPC", "STR")

multiExpr <- vector(mode = "list", length = nSets)
multiExpr[[1]] <- list(data = as.data.frame(t(ctx)))
rownames(multiExpr[[1]]$data) <- names(ctx)

multiExpr[[2]] <- list(data = as.data.frame(t(mb)))
rownames(multiExpr[[2]]$data) <- names(mb)

multiExpr[[3]] <- list(data = as.data.frame(t(hpc)))
rownames(multiExpr[[3]]$data) <- names(hpc)

multiExpr[[4]] <- list(data = as.data.frame(t(str)))
rownames(multiExpr[[4]]$data) <- names(str)

exprSize <- checkSets(multiExpr)
exprSize

traits <- read.csv("traits.csv")

ctx_traits <- brain_key %>%
  filter(tissue == "cortex") %>%
  dplyr::select(c(id, FileID)) %>%
  left_join(traits, by = "id") %>%
  dplyr::select(-c(id, group)) %>%
  mutate(sex = str_replace_all(sex, "F", "1") %>%
           str_replace_all("M", "0")) %>%
  mutate_at(vars(!FileID),
            as.numeric) %>%
  distinct(FileID, .keep_all = T)

mb_traits <- brain_key %>%
  filter(tissue == "midbrain") %>%
  dplyr::select(c(id, FileID)) %>%
  left_join(traits, by = "id") %>%
  dplyr::select(-c(id, group)) %>%
  mutate(sex = str_replace_all(sex, "F", "1") %>%
           str_replace_all("M", "0")) %>%
  mutate_at(vars(!FileID),
            as.numeric) %>%
  distinct(FileID, .keep_all = T)

hpc_traits <- brain_key %>%
  filter(tissue == "hippocampus") %>%
  dplyr::select(c(id, FileID)) %>%
  left_join(traits, by = "id") %>%
  dplyr::select(-c(id, group)) %>%
  mutate(sex = str_replace_all(sex, "F", "1") %>%
           str_replace_all("M", "0")) %>%
  mutate_at(vars(!FileID),
            as.numeric) %>%
  distinct(FileID, .keep_all = T)

str_traits <- brain_key %>%
  filter(tissue == "striatum") %>%
  dplyr::select(c(id, FileID)) %>%
  left_join(traits, by = "id") %>%
  dplyr::select(-c(id, group)) %>%
  mutate(sex = str_replace_all(sex, "F", "1") %>%
           str_replace_all("M", "0")) %>%
  mutate_at(vars(!FileID),
            as.numeric) %>%
  distinct(FileID, .keep_all = T)

allTraits <- rbind(ctx_traits, mb_traits, hpc_traits, str_traits)

multiTraits <- vector(mode = "list", length = nSets)
for (set in 1:nSets){
  setSamples <- rownames(multiExpr[[set]]$data)
  traitRows <- match(setSamples, allTraits$FileID)
  multiTraits[[set]] <- list(data = allTraits[traitRows, -1])
  rownames(multiTraits[[set]]$data) <- setSamples
}

collectGarbage()

nGenes <- exprSize$nGenes
nSamples <- exprSize$nSamples

rownames(multiTraits[[1]]$data) == rownames(multiExpr[[1]]$data)

rownames(multiTraits[[2]]$data) == rownames(multiExpr[[2]]$data)

rownames(multiTraits[[3]]$data) == rownames(multiExpr[[3]]$data)

rownames(multiTraits[[4]]$data) == rownames(multiExpr[[4]]$data)
#Each of the above should return TRUE, indicating the multiExpr and multiTraits \
#lists are identical 

save(multiExpr, multiTraits, nGenes, nSamples, setLabels, shortLabels, exprSize, 
     file = "brain_RNAseq/consensus_WGCNA_input.RData")
