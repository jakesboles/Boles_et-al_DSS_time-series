## GSPA
In between scripts 04 and 05, a gene-set proximity analysis ([Cousins et al., 2023, *Bioinformatics*](https://academic.oup.com/bioinformatics/article/39/1/btac735/6832036?login=true)) was performed on the set of ranked gene lists generated in 04.
The output from this analysis was used in 05. 

Please refer to Henry Cousins's [GSPA](https://github.com/henrycousins/gspa/tree/main) repository for installation and environment setup. I cloned the repository into our lab's partition of the University of Florida' HiPerGator 3.0 so I could run this analysis on a cluster instead of a local machine due to the amount of time I anticipated it would take. 

Before running 05, grab the `.gmt` file containing the gene sets of interest. I used the gene symbols `.gmt` file from [MSigDB's human hallmark gene sets (v.2023.1.Hs)](https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp). Deposit this in `gspa/gene_sets` and name the file however you want. I called it `msigdb_hallmark.gtm` as I was testing other collections. 

In script 04, a set of `.rnk` files are created. I put these files in their own `rnk_files` sub-directory. Then, Python v.3.10 was loaded in the HPC and the GSPA was run with the following commands in a Slurm script:
```
module load python/3.10

cd /PATH/TO/gspa

for i in $(find rnk_files -maxdepth 1 -type f | cut -d '/' -f2 | cut -d '.' -f1); do

INPUT=${i}.rnk
OUTPUT=${i}

python gspa.py \
--rnk_file rnk_files/$INPUT \
--gmt_file gene_sets/msigdb_hallmark.gmt \
--output_folder gspa_output/ \
--results_file $OUTPUT

done
```
After the analysis is finished, move the entire `gspa_output` folder into the `brain_RNAseq` directory of this GitHub repository. 