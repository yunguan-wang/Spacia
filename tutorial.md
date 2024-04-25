## Introduction
This is a simple tutorial for running Spacia on the test dataset manually in linux. 

### Set up
We assume Spacia was already installed. In most cases the python environment is managed by conda or mamba, and R is available through HPC. 
Let's say our python env's name is spacia, and the R version we want to use is 4.1.1.
Load the spacia env.

```
conda activate spacia
module load R/4.1.1
```

### Input data
First, navigate to the `Spacia`'s test data directory.
```
cd [path/to/Spacia/test/input]
```
Two files are essential for running Spacia, one is for gene expression and the other is metadata for each cell.
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

|X|Y|cell_type
| --- | ---| --- |
cell_0|0|0|A
cell_1|0|1|B
cell_2|0|2|B
cell_3|0|3|A

The cell metadata contains spatial coordinates of each cell, as well as its cell type assignment. These three columns and their column names are manditory.

### Test case
In this example, we would like to run Spacia to evaluate how `gene2` and `gene3` in cell_type `B` affects the expression of `gene1` in cell_type `A`.
This can be done in executing the following command.
```
python ../../spacia.py counts.txt spacia_metadata.txt -rc A -sc B -rf gene1 -sf gene2,gene3 -d 5 -m 2000,1000,10,1 -nc 20 -o single_gene_simple_agg
```
Here, `-rc` and `-sc` specified the receiver and sender cell types. 

`-rf gene1` tells Spacia the response gene is `gene1` and `-sf gene2,gene3` indicates the signal genes are gene2 and gene3.

`-d 5` tells Spacia to search for cell type A's neighboring cells within a radius of 5.

`-m` specify the MCMC simulation parameters. In the test example, we set these numbers to be very small to ensure the calculation can finish in minutes. We recommend users to use default parameters for real datasets.

`-nc 20` lets Spacia to use the mean expression of top 20 correlated genes of `gene1` as the response expression. 

Spacia will evaluate the effect of `gene2` on `gene1` and `gene3` on `gene1` simultaneously, and calculate beta values for `gene2` and `gene3`. On the other hand, since there is only one response, only one b value will be calculated, which is for `gene1`, evaluating the dependancy of interactions with `gene1` on cell-cell proximity.

The command should finish in a few minutes, and results and logging information will be saved in the `single_gene_simple_agg` folder. Let's take a look.
```
cd single_gene_simple_agg
ls
```

```
B_and_FDR.csv
gene1_correlated_genes
Interactions.csv
Pathway_betas.csv
model_input
spacia_log.txt
spacia_r.log
```

The output folders contains the summarized results(B_and_FDR.csv, Interactions.csv, Pathway_betas.csv), intermediate results for each response gene(gene1_correlated_genes), and logging information(spacia_log.txt, spacia_r.log).

`B_and_FDR.csv` contains the **b** values of each response gene/pathway (first column) and the associated significance information.

`Pathway_betas.csv` contains the **beta** values representing the interaction between each response gene/pathway (first column) and signal gene/pathway (second columns).

`Interactions.csv` contains the primary instance scores of all receivers in each receiver-sender cell pair (second and third column) for each response-signal interaction (first column). 