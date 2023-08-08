dir1 <- c("brain_RNAseq/", "colon_RNAseq/")
dir2 <- c("csv_outputs", "data_objects", "plots")

for (i in dir1){
  for (j in dir2){
    dir.create(paste0(i, j))
  }
}

dir.create(paste0(dir1[1], "gspa_input"))
dir.create(paste0(dir1[1], "gspa_output"))

renv::restore()