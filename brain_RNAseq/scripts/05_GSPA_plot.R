#This script concatenates the many output files we asked for with the GSPA, tidies the data, and creates a huge plot displaying the results.

library(tidyverse)
library(ggplot2)
library(rlist)
library(data.table)
library(scales)
library(paletteer)
library(stringr)
library(clusterProfiler)

read_plus <- function(flnm) {
  read_csv(flnm) %>%
    mutate(filename = flnm)
}

tbl <- 
  list.files(path = "brain_RNAseq/gspa_output/",
             pattern = "*.csv",
             full.names = T) %>%
  map_df(~read_plus(.))

bonf <- 24

file <- tbl$filename
file <- unlist(strsplit(file, "/", fixed = T)) 
file <- file[c(F, F, T)]
file <- str_remove_all(file, ".csv")
file <- substr(file, 1, 5)

tbl$filename <- file 

tbl <- tbl %>%
  mutate(FDR = FDR * bonf) %>%
  separate(filename, into = c("tissue", "group"), sep = "_") %>%
  mutate(tissue = factor(tissue,
                         levels = c("c", "h", "m", "s"),
                         labels = c("Cortex", "Hippocampus", 
                                    "Midbrain", "Striatum")),
         group = str_remove_all(group, "gr") %>%
           factor(levels = 1:6,
                  labels = c("5d DSS", "7d DSS", "7d DSS + 2d H2O",
                             "7d DSS + 5d H2O", "7d DSS + 7d H2O",
                             "7d DSS + 14d H2O")),
         `Gene Set` = str_replace_all(`Gene Set`, "_", " ") %>%
           str_remove_all("HALLMARK "))

tbl %>%
  na.omit() %>%
  filter(FDR < 0.05) %>%
  mutate(reg = case_when(NES > 0 ~ "Up-regulated",
                         NES < 0 ~ "Down-regulated")) %>%
  ggplot(aes(x = group, y = factor(`Gene Set`,
                                   levels = rev(levels(factor(`Gene Set`)))))) +
  facet_grid(cols = vars(tissue), rows = vars(factor(reg,
                                                     levels = c("Up-regulated",
                                                                "Down-regulated"))),
             scales = "free_y", space = "free_y") + 
  geom_point(aes(size = abs(NES), color = FDR)) +
  scale_color_paletteer_c(palette = "grDevices::ag_Sunset",
                          direction = 1) +
  scale_x_discrete(labels = label_wrap(9)) +
  labs(size = "Normalized\nenrichment score\n(absolute value)",
       color = "Bonferroni-\nadjusted FDR") +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(color = "white", face = "bold", size = 16),
    axis.title = element_blank(),
    axis.text = element_text(color = "black", face = "bold")
  )

ggsave("brain_RNAseq/plots/gspa_plot.png",
       units = "in", dpi = 600,
       height = 10, width = 20)
