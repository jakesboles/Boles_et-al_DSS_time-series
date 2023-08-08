# Boles_et-al_DSS_time-series

**This repository is currently under construction**

### This repository contains code to reproduce analysis results and visualizations for data in:
### **A leaky gut dysregulates gene networks in the brain associated with immune activation, oxidative stress, and myelination in a mouse model of colitis (*bioRxiv*)**
Jake Sondag Boles<sup>1</sup>, Maeve E. Krueger, Janna E. Jernigan, Cassandra L. Cole, Noelle K. Neighbarger, Oihane Uriarte Huarte, & Malú Gámez Tansey.

<sup><sup>1</sup> Analysis lead and contact (jake.boles@ufl.edu)</sup>

## Data:
### Original bulk RNA-sequencing data:
This can be accessed via the NCBI GEO ([GSE239820](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239820)). We have included:
1. Raw FASTQ files from paired-end reads (two files per sample) 
2. Raw transcript counts after read mapping and counting. These data are in "raw_counts.csv". The key can be found in the metadata at the GEO entry page. They will also be scraped from the link in `01_counts_cleaning.R` for both colon and brain pipelines. 
3. Processed and normalized counts after the removal of low-quality genes and samples. These are the data used to arrive at the conclusions in the paper above. These data are separated into colon ("colon_clean_VST_counts.csv") and brain ("brain_clean_VST_counts.csv"). The sample key can be found in the metadata at the GEO entry page. This repository contains the code that was used to generate these cleaned counts from the raw counts. 

### Publicly available datasets used here:
| Dataset | Publication | DSS schedule | Colon regions examined |
| :-----: | :---------: | :----------: | :--------------------: |
| [GSE131032](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131032) | [Czarnewski et al., 2019 <br> *Nature Communcations*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6598981/) | 2.5% DSS for 7d, followed by <br> 7d of recovery. Tissue taken <br> at days 0, 2, 4, 6, 8, 10, 12, 14 | Mid-colon |
| [GSE168053](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168053)<sup>#</sup> | [Liu et al., 2022 <br> *Gastroenterology*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9402284/) | 3% DSS for 6d, followed by <br> 14d of recovery. Tissue taken <br> at days 0, 3, 6, 9, 14, and 20 | Mid-colon ("proximal"), <br> distal, anal, <br> squamous neo-epithelium |
| [GSE210405](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210405) | NA | 3% DSS for 6d, with tissue <br> harvested at days 0, 2, 4, and 6 | Mid-colon | 

<sup><sup>#</sup>Only mid-colon (termed "proximal") and distal colon samples were used in our study to compare with our distal and proximal samples.

### Other original data:
This includes data procured from lower-throughput assays used in the paper, including flow cytometry, immunoblotting, RT-qPCR, ELISAs, MesoScale Discovery assays, and daily disease monitoring.

## Code: 
This project was run in R v4.2. It uses the `renv` package for package management.
All packages were loaded and up to date on July 31, 2023. The library was constructed to be compatible with R v.4.2, using Bioconductor v.3.16. 

**NOTE: If you intend to run this pipeline yourself, it may be prudent to download the repository to your institution's HPC if available.** A few scripts including the construction of the topological overlap matrix for co-expression module identification and module preservation meta-analyses were done on the University of Florida's HiPerGator 3.0, with 85GB RAM and 11 CPUs. These steps are very computationally intensive, although I did not evaluate if they would run to completion on a less powerful machine. 

I might also gingerly recommend *avoiding* a Mac OS. I have had trouble loading dependencies for key packages such as `DESeq2` through `renv`, and it seems to be a problem with the Fortran compiler, which may be very difficult to address especially if using a university-owned Mac machine. Making this pipeline accessible to all systems is currently a work in progress.  

To get started, open a new R session and run:
```
if (!require(usethis)) {
install.packages("usethis")
}

library(usethis)

usethis::use_course(
  'jakesboles/Boles_et-al_DSS_time-series',
  destdir = 'PATH/TO/DIRECTORY')
```
This uses the `usethis` package, and it will be installed on your machine if it is not found already. Specify the destination directory appropriately. A new folder in the destination directory will be created and the new RStudio session should launch.  

Then, initialize the project by running `00_initialize_project.R` in the newly created R Project. This will create sub-folders to store output from the pipeline and install all version-locked packages. You may be prompted to confirm package installation in the console. This may take some time! 

