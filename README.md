# Boles_et-al_DSS_time-series
### Code to reproduce analysis results and visualizations for data in:
### **An unbiased examination of a leaky gut's impact on the brain in a mouse model of colitis (*bioRxiv*)**
Boles, J. S., Krueger, M. E., Jernigan, J. E., Cole, C. L., Neighbarger, N. K., Uriarte Huarte, O., & Tansey, M. G.

##Data:
Bulk RNA-sequencing data can be accessed via the NCBI GEO (insert GSE number when submission is completed). We have included:
1. Raw FASTQ files from paired-end reads (two files per sample) 
2. Raw transcript counts after read mapping and counting. These data are in "raw_counts.csv"(insert link). The key can be found in the metadata at the GEO entry page. 
3. Processed and normalized counts after the removal of low-quality genes and samples. These are the data used to arrive at the conclusions in the paper above. These data are separated into colon ("colon_VST_counts_3.csv"(insert link)) and brain ("brain_VST_counts_3.csv"(insert link)). The sample key can be found in the metadata at the GEO entry page. 
