# Boles_et-al_DSS_time-series

**This repository is currently under construction**

### This repository contains code to reproduce analysis results and visualizations for data in:
### **An unbiased examination of a leaky gut's impact on the brain in a mouse model of colitis (*bioRxiv*)**
Boles, J. S.<sup>1</sup>, Krueger, M. E., Jernigan, J. E., Cole, C. L., Neighbarger, N. K., Uriarte Huarte, O., & Tansey, M. G.

<sup><sup>1</sup> Analysis lead and contact (jake.boles@ufl.edu)</sup>

## Data:
### Original bulk RNA-sequencing data:
This can be accessed via the NCBI GEO (insert GSE number when submission is completed(insert link)). We have included:
1. Raw FASTQ files from paired-end reads (two files per sample) 
2. Raw transcript counts after read mapping and counting. These data are in "raw_counts.csv"(insert link). The key can be found in the metadata at the GEO entry page. 
3. Processed and normalized counts after the removal of low-quality genes and samples. These are the data used to arrive at the conclusions in the paper above. These data are separated into colon ("colon_VST_counts_3.csv"(insert link)) and brain ("brain_VST_counts_3.csv"(insert link)). The sample key can be found in the metadata at the GEO entry page. 

### Publicly available datasets used here:
| Dataset | Publication | DSS schedule | Colon regions examined |
| :-----: | :---------: | :----------: | :--------------------: |
| [GSE131032](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131032) | [Czarnewski et al., 2019 <br> *Nature Communcations*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6598981/) | 2.5% DSS for 7d, followed by <br> 7d of recovery. Tissue taken <br> at days 0, 2, 4, 6, 8, 10, 12, 14 | Mid-colon |
| [GSE168053](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168053)<sup>a</sup> | [Liu et al., 2022 <br> *Gastroenterology*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9402284/) | 3% DSS for 6d, followed by <br> 14d of recovery. Tissue taken <br> at days 0, 3, 6, 9, 14, and 20 | Mid-colon ("proximal"), <br> distal, anal, <br> squamous neo-epithelium |
| [GSE210405](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210405) | NA | 3% DSS for 6d, with tissue <br> harvested at days 0, 2, 4, and 6 | Mid-colon | 

<sup><sup>a</sup>Only mid-colon (termed "proximal") and distal colon samples were used in our study to compare with our distal and proximal samples.

### Other original data:
This includes data procured from lower-throughput assays used in the paper, including flow cytometry, immunoblotting, RT-qPCR, ELISAs, MesoScale Discovery assays, and daily disease monitoring.

## Code: 
This project was run in R v4.2-4.3. It uses the `renv` package for package management.
All packages were loaded and up to date on July 31, 2023. The library was constructed to be compatible with R v.4.2, using Bioconductor v.3.16. 
To clone this Github repository and adopt the version-controlled library, start a new RStudio session and run the following:
```
library(usethis)

create_from_github(
  'jakesboles/Boles_et-al_DSS_time-series',
  destdir = 'PATH/TO/DIRECTORY',
  fork = F,
  protocol = "https"
)
```
Note that this requires the `usethis` package, which you may need to install with `install.packages("usethis")` before executing the above code if you don't have it already. 
This command will open a new RStudio Project file in the directory you specified. To install all packages:
```
renv::restore()
```
When prompted, input `y` in the console to accept changes to the renv.lock file. This could take some time if you are initializing the project for the first time!

