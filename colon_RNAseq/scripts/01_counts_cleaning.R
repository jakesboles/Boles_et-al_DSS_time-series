library(DESeq2)
library(tidyverse)

#Grab raw counts file from GEO's https link and convert to data frame

con <- gzcon(url("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE239820&format=file&file=GSE239820%5Fraw%5Fcounts%2Ecsv%2Egz"))
csv <- readLines(con)
counts <- read.csv(textConnection(csv))