## Introduction
This is a tutorial for manually running the Spacia R interface using the test dataset in linux. 
The R interface runs Spacia in a quality-optimized mode that does not support some of the options in the python interface.

### Set up
The R interface does not require python, but it does require the same core R packages as well as some additional packages as detailed in the README.
If only runnning the R interface, then install/load the correct version of R (and associated libraries/compilers if needed) and install the packages with:

```
#example on a generic HPC environment
module load R/4.1.1 gcc/12.2.0

#start R session
R

#install core packages
install.packages(c('coda', 'ggmcmc', 'Rcpp', 'RcppArmadillo', 'rjson'))

#install R interface specific packages
install.packages(c('optparse', 'filelock', 'ggplot2', 'patchwork', 'scales', 'gridExtra', 'dplyr'))

#exit after all packages are installed
Ctrl+D
```

### Input data
The R interface can use the same input data as the python interface. Exit R and navigate to `Spacia`'s test data directory.
```
cd [path/to/Spacia/test/input]
```
Two files are essential for running Spacia, one is for gene expression and the other is the location metadata for each cell.
Let's take a look at the expression data.
```
head -n 5 counts.txt | cut -d $'\t' -f 1-6
```

| | gene1 | gene2 | gene3 | gene4 | gene5 |
| --- | ---| --- | ---| --- | ---|
| cell_0 | 1.058939863 | 2.142738423 | 1.357627858 | 0.942371598 | 1.523920809 |
| cell_1 | 0.973380302 | 2.423213296 | 1.428899385 | 1.208159298 | 1.165738144 |
| cell_2 | 0.75562576 | 2.164291889 | 1.067379324 | 1.456489891 | 1.466352688 |
| cell_3 | 0.820909431 | 2.011477308 | 1.249298614 | 1.181445221 | 2.132801071 |

As we can see, the gene expression data is a cell-by-gene matrix, where the first column is cell names, and the first row is the gene names.
This test data is a matrix of 2,844 cells by 100 simulated genes.

Then, let's take a look at the metadata.
```
head -n 5 spacia_metadata.txt
```

| |X|Y|cell_type
| --- | ---| --- | ---|
cell_0|0|0|A
cell_1|0|1|B
cell_2|0|2|B
cell_3|0|3|A

The cell metadata contains spatial coordinates of each cell, as well as its cell type assignment. These three columns and their column names are manditory. 
Cell names must be consistent between the two input files. Also make sure that the coordinates of the cells are consistent and in the correct units; 
we highly recommend using microns as the unit of measure for distances.

### Test case
In this example, we would like to run Spacia to evaluate how genes in cell_type `A` affects the expression of `gene2` in cell_type `B`. 

First, we need to find the necessay cutoff values. Unlike in the python interface, these values need to be speficied manually in the current version.
Simply run Spacia without specifying any of the required cutoffs and plots will be generated to aid in determining the values.
```
Rscript spacia.R \
	-x counts.txt \
	-m spacia_metadata.txt \
	-a spacia \
	-r B -s A -g gene2 \
	-d 2 \
	-o r_test/
```
Here, `-r` and `-s` specifies the receiver and sender cell types. 

For sending genes, the R interface always aggregates all genes using PCA, 
which is the same behavior as passing `pca` to `--sender_features` in the python interface. 

`-g` specifies the receiving gene. There is no option to aggregate the receiving genes.
If mulitple receiving genes are supplied by passing a file, 
then each gene will be processed in sequence in seperate Spacia runs.

`-d 2` tells Spacia to consider "A" cells within 2 distance radius of "B" cells as potential interacting cells. 
This small value is intended to match the synthetic test data, and so should not be used for real data.

`-o` specifies the output prefix. Inputs ending with "/" are intereted as directories, 
therefore `r_test` will be created to contain the outputs.

Spacia will generate a cache file `A-B_cache.RData` and a plots pdf `gene2_cutoffPlot.pdf` under `r_test`.
Open `gene2_cutoffPlot.pdf` and determine suitable cutoffs according to the instructions in the README. 
For this case, we can use 0.252 as the receiving gene correlation cutoff and 0.091 as the receiving gene quantile cutoff.

We can then run Spacia using the cutoffs:
```
Rscript spacia.R \
	-x counts.txt \
	-m spacia_metadata.txt \
	-a spacia \
	-r B -s A -g gene2 \
  -q 0.091 -u 0.252 \
	-d 2 -l 5000 -w 2500 \
	-o r_test/
```

`-q` and `-u` are the receiving gene correlation and quantile cutoff values.

`-l` and `-w` specify the MCMC simulation parameters. 
In the test example, we set these numbers to be very small to ensure the calculation can finish in minutes. 
We recommend users to use default parameters for real datasets.

Alternatively, we can organize the cutoffs into a table (csv or tsv):
```
head gene_cutoffs_A-B.txt
```
|gene_name|cor_cutoff|quantile_cutoff
| --- | ---| --- |
gene2|0.252|0.091

The gene cutoffs file can then be used instead of manually inputing the cutoffs:
```
Rscript spacia.R \
	-x counts.txt \
	-m spacia_metadata.txt \
	-a spacia \
	-r B -s A \
  -t gene_cutoffs_A-B.txt \
	-d 2 -l 5000 -w 2500 \
	-o r_test/
```

In this case, we can omit `-g`, which causes Spacia to attempt to sequentially process all genes listed in `gene_cutoffs_A-B.txt`, 
which is only `gene2` in this example. This can be useful in real analyses with many receiving genes. 
If using the same storage system and output directory, multiple Spacia instances can be run simultaneously using the same command. 
Otherwise, split the table for `-t` accordingly. Users can also specify a list of receiving genes by passing a text file to `-g`.

For this example Spacia will evaluate the effect all genes from `A` cells to `gene2` in `B` cells simultaneously, 
and calculate beta values for each gene using its PCA weights. 
On the other hand, since there is only one response, only one b value will be calculated, which is for `gene2`, 
evaluating the dependancy of interactions with `gene1` on cell-cell proximity.

The command should finish in a few minutes, and results and logging information will be saved in the `r_test` folder. Let's take a look.
```
cd r_test
ls
```

```
A-B_cache.RData
A-B_gene2_pip.csv  
A-B_gene2_betas.csv
A-B_gene2.RData
gene2_cutoffPlot.pdf
```

The output folders contains the summarized results (A-B_gene2_betas.csv, A-B_gene2_pip.csv), 
a combined results data file (A-B_gene2.RData), and the previously generated plots (gene2_cutoffPlot.pdf).

`A-B_gene2_betas.csv` contains both the the **beta** values representing the 
interaction between each signal gene/pathway (first columns) and response gene/pathway (second column),
and the **b** values of each response gene/pathway (second column) along with the associated significance information.
Since there is only one receiving gene per run, the **b** values are the same for all rows. 

`A-B_gene2_pip.csv` contains the primary instance scores of all receivers in each receiver-sender cell pair (first and second column) 
for each response-signal interaction (third column). 
