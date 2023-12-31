library(readxl)
library(tidyverse)

attr <- read.csv("traits.csv")
attr <- attr[,c(1:3)]
attr$group <- factor(attr$group,
                     levels = c("Untreated", "5d_DSS", "7d_DSS", "7d_DSS_2d_H2O",
                                "7d_DSS_5d_H2O", "7d_DSS_7d_H2O", "7d_DSS_14d_H2O"))

key <- read_xlsx("IGNORE/9441.xlsx")
key <- key %>%
  separate(col = "Sample ID",
           into = c("id", "tissue"), sep = "_") %>%
  left_join(attr, by = "id") %>%
  mutate(tissue = case_when(
    tissue == "C" ~ "cortex",
    tissue == "S" ~ "striatum",
    tissue == "M" ~ "midbrain",
    tissue == "H" ~ "hippocampus",
    tissue == "D" ~ "distal_colon",
    tissue == "P" ~ "proximal_colon"
  )) %>%
  distinct(FileID, .keep_all = T)

gut <- key[c(266:nrow(key)), ]
brain <- key[c(1:265), ]

write.csv(brain, file = "brain_RNAseq/RNAseq_brain_metadata.csv",
          row.names = F)
write.csv(gut, file = "colon_RNAseq/RNAseq_colon_metadata.csv",
          row.names = F)
