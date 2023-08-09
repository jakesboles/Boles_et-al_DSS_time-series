#This script creates the input files for the consensus WGCNA
#Note that we skip over a lot of the analysis performed in brain, including DE analysis and GSPA
#This is simply because those analyses are not novel in colon, and we wanted to leverage this fact by using similar colon datasets

library(WGCNA)
library(tidyverse)
library(DESeq2)

load("colon_RNAseq/data_objects/deseq2_objects.RData")

dc <- vsd[, colnames(vsd) %in% gut_key[grep("distal_colon", gut_key$tissue), ]$FileID]
pc <- vsd[, colnames(vsd) %in% gut_key[grep("proximal_colon", gut_key$tissue), ]$FileID]

nSets <- 2
setLabels <- c("Distal colon", "Proximal colon")
shortLabels <- c("Distal", "Proximal")

multiExpr <- vector(mode = "list", length = nSets)
multiExpr[[1]] <- list(data = as.data.frame(t(dc)))
names(multiExpr[[1]]$data) <- names(dc)

multiExpr[[2]] <- list(data = as.data.frame(t(pc)))
names(multiExpr[[2]]$data) <- names(pc)

exprSize <- checkSets(multiExpr)
exprSize

traits <- read.csv("traits.csv")

dc_traits <- gut_key %>%
  filter(tissue == "distal_colon") %>%
  dplyr::select(c(id, FileID)) %>%
  left_join(traits, by = "id") %>%
  dplyr::select(-c(id, group)) %>%
  mutate(sex = str_replace_all(sex, "F", "1") %>%
           str_replace_all("M", "0")) %>%
  mutate_at(vars(!FileID),
            as.numeric) %>%
  distinct(FileID, .keep_all = T)

pc_traits <- gut_key %>%
  filter(tissue == "proximal_colon") %>%
  dplyr::select(c(id, FileID)) %>%
  left_join(traits, by = "id") %>%
  dplyr::select(-c(id, group)) %>%
  mutate(sex = str_replace_all(sex, "F", "1") %>%
           str_replace_all("M", "0")) %>%
  mutate_at(vars(!FileID),
            as.numeric) %>%
  distinct(FileID, .keep_all = T)

allTraits <- rbind(dc_traits, pc_traits)

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
#Each of the above should return TRUE, indicating the multiExpr and multiTraits \
#lists are identical 

save(multiExpr, multiTraits, nGenes, nSamples, setLabels, shortLabels, exprSize, 
     file = "colon_RNAseq/consensus_WGCNA_input.RData")