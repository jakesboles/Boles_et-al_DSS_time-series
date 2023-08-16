#This script annotates all colon consensus modules with anRichment 
#These results are used to identify biologically interesting modules
#Horizontal bar charts will also be created to display the results

library(WGCNA)
library(anRichment)
library(ggplot2)
library(tidyverse)
library(stringr)
library(RColorBrewer)
library(paletteer)
library(forcats)
library(scales)
library(cowplot)
library(ggpubr)

options(stringsAsFactors = F)

load("colon_RNAseq/consensus_WGCNA_input.RData")
load("colon_RNAseq/data_objects/consensus_network.RData")

#Enrichment analysis------------------------------------------------------------
symbol <- colnames(multiExpr[[1]]$data)

entrez <- convert2entrez(organism = "Mus musculus", symbol = symbol)

GOcollection <- buildGOcollection(organism = "mouse")
knownGroups(GOcollection)
GO.BPcollection <- subsetCollection(GOcollection, tags = "GO.BP")

#subsetting GO.BP collection to remove huge gene lists----
idx <- vector(mode = "logical", length = length(GO.BPcollection$dataSets))

max_genes <- 500

for (i in seq_along(GO.BPcollection$dataSets)){
  idx[i] <- nrow(GO.BPcollection$dataSets[[i]]$data) < max_genes
}

small_go <- GO.BPcollection$dataSets[idx]
{
  dataSets <- list()
  
  for (i in seq_along(small_go)){
    go_entrez <- small_go[[i]]$data$Entrez
    go_evidence <- small_go[[i]]$data$evidence
    go_source <- small_go[[i]]$data$source
    
    dataSets[i] <- newGeneSet(
      geneEntrez = go_entrez,
      geneEvidence = go_evidence,
      geneSource = go_source,
      ID = small_go[[i]]$ID,
      name = small_go[[i]]$name,
      description = small_go[[i]]$description,
      source = small_go[[i]]$source,
      organism = "mouse",
      internalClassification = small_go[[i]]$internalClassification,
      groups = small_go[[i]]$groups,
      lastModified = small_go[[i]]$lastModified
    )
  }
  
  names(dataSets) <- names(small_go)
}

small_go_group <- newGroup(name = "GO.BP",
                           description = paste0("GO-BP knowledgebase gene sets with fewer than ", max_genes, " genes"),
                           source = "GO")
{
  smallGOcollection <- newCollection(groups = list(small_go_group))
  smallGOcollection$dataSets <- vector(mode = "list", length = length(dataSets))
  
  for (i in seq_along(dataSets)){
    smallGOcollection$dataSets[[i]]$data <- data.frame(Entrez = dataSets[[i]]$Entrez,
                                                       evidence = dataSets[[i]]$evidence,
                                                       source = dataSets[[i]]$source)
    
    names(smallGOcollection$dataSets)[[i]] <- names(dataSets)[[i]]
    
    smallGOcollection$dataSets[[i]]$ID <- small_go[[i]]$ID
    
    smallGOcollection$dataSets[[i]]$name <- small_go[[i]]$name
    
    smallGOcollection$dataSets[[i]]$shortName <- small_go[[i]]$shortName
    
    smallGOcollection$dataSets[[i]]$description <- small_go[[i]]$description
    
    smallGOcollection$dataSets[[i]]$source <- small_go[[i]]$source
    
    smallGOcollection$dataSets[[i]]$organism <- small_go[[i]]$organism
    
    smallGOcollection$dataSets[[i]]$internalClassification <- small_go[[i]]$internalClassification
    
    smallGOcollection$dataSets[[i]]$groups <- small_go[[i]]$groups
    
    smallGOcollection$dataSets[[i]]$lastModified <- small_go[[i]]$lastModified
    
    smallGOcollection$dataSets[[i]]$alternateNames <- small_go[[i]]$alternateNames
    
    smallGOcollection$dataSets[[i]]$externalDB <- small_go[[i]]$externalDB
    
    smallGOcollection$dataSets[[i]]$externalAccession <- small_go[[i]]$externalAccession
    
    smallGOcollection$dataSets[[i]]$webLink <- small_go[[i]]$webLink
    
    smallGOcollection$dataSets[[i]]$weightIndex <- small_go[[i]]$weightIndex
    
    smallGOcollection$dataSets[[i]]$type <- small_go[[i]]$type
  }
}

#Build other collections
biosysCollection <- BioSystemsCollection("mouse")
knownGroups(biosysCollection)
keggCollection <- subsetCollection(biosysCollection, tags = "KEGG",
                                   matchComponents = "groups")
biocycCollection <- subsetCollection(biosysCollection, tags = "BIOCYC",
                                     matchComponents = "groups")
reactomeCollection <- subsetCollection(biosysCollection, tags = "REACTOME",
                                       matchComponents = "groups")

msdbCollection <- MSigDBCollection(file = "msigdb_v2023.1.Mm.xml/msigdb_v2023.1.Mm.xml", organism = "mouse")
mhCollectrion <- subsetCollection(msdbCollection,
                                  tags = "MSigDB MH: NA")
wikiCollection <- subsetCollection(msdbCollection,
                                   tags = "MSigDB M2: NA - CP:WIKIPATHWAYS")
biocartaCollection <- subsetCollection(msdbCollection,
                                       tags = "MSigDB M2: NA - CP:BIOCARTA")

combinedCollection <- mergeCollections(smallGOcollection, biocartaCollection,
                                       keggCollection, mhCollectrion,
                                       reactomeCollection, wikiCollection)

combined <- enrichmentAnalysis(classLabels = moduleColors,
                               identifiers = entrez,
                               refCollection = combinedCollection,
                               useBackground = "given",
                               threshold = 0.05,
                               thresholdType = "FDR",
                               getOverlapSymbols = T,
                               maxReportedOverlapGenes = 300)

write.csv(combined$enrichmentTable, 
          file = "colon_RNAseq/csv_outputs/anRichment_output.csv",
          row.names = F)

df <- combined$enrichmentTable

#Preparing data frame for visualizations----------------------------------------
extract_transform_value <- function(string) {
  # Check if "WIKIPATHWAYS" is present
  if (grepl("WIKIPATHWAYS", string)) {
    transformed_value <- gsub(".*WIKIPATHWAYS.*", "WIKIPATHWAYS", string)
  } else {
    transformed_value <- string
  }
  
  transformed_value
}
# Apply the function to the strings column
df$inGroups <- sapply(df$inGroups, extract_transform_value)

extract_transform_value <- function(string) {
  # Check if "MH" is present
  if (grepl("MH", string)) {
    transformed_value <- gsub(".*MH.*", "MSIGDB", string)
  } else {
    transformed_value <- string
  }
  
  transformed_value
}
# Apply the function to the strings column
df$inGroups <- sapply(df$inGroups, extract_transform_value)

extract_transform_value <- function(string) {
  # Check if "CGP" is present
  if (grepl("CGP", string)) {
    transformed_value <- gsub(".*CGP.*", "MSIGDB CGP", string)
  } else {
    transformed_value <- string
  }
  
  transformed_value
}
# Apply the function to the strings column
df$inGroups <- sapply(df$inGroups, extract_transform_value)

extract_transform_value <- function(string) {
  # Check if "MH" is present
  if (grepl("BIOCARTA", string)) {
    transformed_value <- gsub(".*BIOCARTA.*", "BIOCARTA", string)
  } else {
    transformed_value <- string
  }
  
  transformed_value
}
# Apply the function to the strings column
df$inGroups <- sapply(df$inGroups, extract_transform_value)

extract_transform_value <- function(string) {
  # Check if "GO.BP" is present
  if (grepl("", string)) {
    transformed_value <- gsub(".*GO.BP.*", "GO.BP", string)
  } else {
    transformed_value <- string
  }
  
  transformed_value
}
# Apply the function to the strings column
df$inGroups <- sapply(df$inGroups, extract_transform_value)

df$dataSetName <- gsub("\\s*\\([^)]+\\)", "", df$dataSetName)

convert_pathway_names <- function(string) {
  # Split the string by '|' or underscores
  words <- unlist(strsplit(string, "[|_]"))
  
  # Convert each word to uppercase
  words_upper <- toupper(words)
  
  # Join the words with spaces in between
  converted_string <- paste(words_upper, collapse = " ")
  
  converted_string
}
df$dataSetName <- sapply(df$dataSetName, convert_pathway_names)

df$geneSet <- paste(df$inGroups,
                    df$dataSetName,
                    sep = ": ")
unique(df$inGroups)
#Plotting-----------------------------------------------------------------------
lollipop <- function(module, titlecolor = NULL, wrap = NULL){
  d <- data %>%
    mutate(fracOfEffectiveClassSize = 100 * fracOfEffectiveClassSize) %>%
    subset(class == module) %>%
    distinct(geneSet, .keep_all = T) %>%
    #mutate(dataSetName = paste(inGroups,
    #                           dataSetName,
    #                           sep = ": ")) %>%
    arrange(pValue) %>%
    mutate(geneSet = factor(geneSet,
                            levels = geneSet)) %>%
    slice_head(n = 10)
  
  title <- str_to_title(module)
  
  if(is.null(titlecolor)) titlecolor <- module
  if(is.null(wrap)) wrap <- 20
  
  d %>%
    ggplot(aes(x = -log10(FDR), y = geneSet)) + 
    geom_col(fill = module, color = "black") +
    #geom_point(aes(color = nCommonGenes),
    #           size = 12) + 
    scale_y_discrete(labels = label_wrap(wrap),
                     limit = rev) +
    scale_color_paletteer_c("grDevices::Spectral",
                            direction = 1) +
    geom_vline(xintercept = -log10(0.05),
               linetype = "twodash",
               linewidth = 1.5) +
    labs(x = "-log10(FDR)",
         color = "# of\noverlapping\ngenes") + 
    ggtitle(paste0(title)) +
    guides(color=guide_colourbar(barheight=20,label.position="right",
                                 barwidth = 2, title.hjust = 0.5,
                                 title.position = "top")) +
    theme_bw(base_size = 12) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 20, color = "black", face = "bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          axis.text.x = element_text(size = 14, color = "black"),
          axis.title.x = element_text(size = 16, color = "black"),
          plot.title = element_text(size = 20, color = "black", face = "bold",
                                    hjust = 0.5))
}

saveplot <- function(filename){
  ggsave(paste0("colon_RNAseq/plots/", filename),
         units = "in", dpi = 600,
         height = 9, width = 8)
}

fig2 <- c("blueviolet", "darkseagreen3", "orange", "darkorange", "red", "firebrick3")
fig2_wrap <- c()
sfig4 <- c("coral2", "coral3", "white", "thistle2")
sfig4_wrap <- c()
sfig5 <- c("brown", "coral", "darkolivegreen4", "lightpink3", "coral1", "indianred4", "darkslateblue")
sfig5_wrap <- c()
sfig6 <- c("floralwhite", "lightcoral", "salmon", "tan")
sfig6_wrap <- c()

figs <- list(fig2, sfig4, sfig5, sfig6)
wraps <- list(fig2_wrap, sfig4_wrap, sfig5_wrap, sfig6_wrap)

for (i in seq_along(figs)){
  for (j in seq_along(figs[[i]])){
    lollipop(figs[[i]][j], wrap = wraps[[i]][j])
    
    saveplot(paste0(deparse_substitute(figs[[i]]), "_", figs[[i]][j], ".png"))
  }
}
