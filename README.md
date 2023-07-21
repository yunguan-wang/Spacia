# Spacia: Inferring Cell-cell Interaction from Spatial Transcriptomics Data


## Introduction
Multicellular organisms heavily rely on cell-cell interactions to effectively coordinate and regulate various biological processes, ensuring the normal functioning of the organism. Spacia models and evaluates cell-cell interactions from spatial transcriptomic data (SRT). This model uses cell-cell proximity as a constraint and prioritizes cell-cell interactions that cause a downstream change in the cells. Spacia employs a Bayesian multi-instance learning (MIL) framework to assess intercellular communication from between cells and their neighbors.

## Graphical abstract
Insert Fig2b high res png file 

## Installation
### Dependency
R: `R>=4.0`, `coda`, `ggmcmc`, `rcpp`, `RcppArmadillo`.

Python: `R>=3.8`, `matplotlib`, `pandas`, `scipy`, `scikit-learn`.
We strongly recommend using conda to manage the installation of all dependencies. To do this, simply run:

```
conda create --name spacia
conda activate spacia
# conda config --add channels conda-forge ##If you do not have this channel added before#
conda install r-base=4.0 r-coda r-ggmcmc r-rcpp r-RcppArmadillo
conda install python=3.8 pandas scikit-learn numpy scipy matplotlib
```
Then, download this repo.
```
git clone [repo_path]
```

The total installation time is around 10 mintunes. If error occuors, please upgrade pip and try again.

#### Note:

## Test installation
Todo

## Usage

### Definition of terms used in Spacia
**Interaction**: Relationship between a pair **interactants** that potentially leads to downstream signalling events in cells. The **interactant** can be a gene or a geneset.

**Signal**: The **interactant** in the **interaction** that is causing downstream signaling events.

**Response**: The **interactant** in the **interaction** whose expression is changed as the result of activities from **Signal**.

**Sender**: A cell where the **Signal** is expressed.

**Receiver**: A cell where the **Response** is expressed. 

**Neighborhood**: A regions centering around each **Receivers** that contains **Senders** of interest.

### Quick start
Once the input data have been processed into the supported format, the full Spacia workflow can be run by calling the `Spacia.py` script. It evaluates interactions within the context of cell neighborhoods, where the ‘**receiver**’ cells are the cells of interest, and the cells from the neighborhood are referred to as "**sender**" cells. The **interactant** expressed in the receiver cells, through which the interactions are to be studied, are referred to as "**Response**", while the **interactant** expressed in the sender cells that potentially influence the responder genes are called signal “**Signal**".

```
python [path/to/spacia_job.py] counts.txt cell_metadata.txt -rc celltype1 sc celltype2 -rf gene1 sf gene2
```

Here, `counts.txt` is a cell-by-gene matrix in TX format. We expect the data to be normliazed, if not, CPM normalization will be used.

`cell_metadata.txt` is a cell_by_metadata matrix in txt format in TXT format. Must contains `X` and `Y` columns for coordinates, and a `cell_type` columns, refering to the group designation of cells, is needed if '-rc' or '-sc' parameter are given.

`-rc` and `-sc` refer to **receiver** cells and **sender** cells, respectively.

`-rf` and `-sf` refer to **Response** and **Signal** features. Here they are in forms of single genes. Spacia can also take pathways in the format of a list of genes as inputting features. 
### List of other important parameters
`--dist_cutoff` or `-d`: The euclidean distance deifining the radius of the neighborhood around each receiver cell.

### Processing **interactant** expression
**Spacia employs several different workflows to calculate **interactant** expression in cells, aiming to handle use cases of disfferent purposes. The behavior is controled largely by the `--receiver_features` and `--sender_features` paramerters, and a few others by a lesser extent. **
* When the **interactant** is a single gene, Spacia can try to mitigate noises associated with gene expression in SRT data by considering the expression of highly correlated genes (by absolute Pearson correlation values). This behavior can be turned off by passing the `--corr_agg` keyword. The number of highly correlated genes to consider can be changed by passing the desired number to the `--num_corr_genes` keyword. The new expression value considering these correlated genes can be calculated as the weighted average of there expression, whereas the weights are the Pearson correlation coefficients with the gene of interest. In cases where only the positively correlated genes should be considered, spacia will only include the top positively correlated genes to calculate the expression of the **interactant**. This behavior can be set by passing `--corr_agg_method simple`.
* When the **interactant** contains serveral genes, Spacia will no longer use correlation based aggregation, instead, the average of the inputs genes will be calculated and used as the expression of the interactant.
The list of genes can be passed as a string separted by "|", e.g., 'CD3E|CD4|CD8A'. It can be also passed by a csv files, with each genelist as a separated row, where the first element is the name of the genelist. These paramerters should be passed to `--receiver_features` or `--sender_features`.
* Spacia can also be run in two unsupervised modes where the **interactant** is not provided. In the first unsupervised mode, spacia will transform the SRT data using the first 20 principal components, and use the transformed dimensions as **interactants**. ****This mode is not recommended for response genes, as the interactions predicted in this way are difficult to interpret.**** This mode can be set by passing `pca` to `--receiver_features` or `--sender_features`. In the second unsupervised mode, spacia will cluster the genes in the SRT data using hierarchical clustering and use the expression value of each cluster centroid as the **interactions**. This mode can be set by not passing any parameters to `--receiver_features` or `--sender_features`.

### How to use a custom list of cells as **receiver** or **sender**

### List of Parameters
```
Todo
```

### For advanced users
For users who want to directly access the core of spacia and perform more flexible analyses, we provide an example R scipt that showcases the few key steps. But please regard the codes in this R script as examples and remember to customize everything according to your needs/datasets. Our analysis codes of the prostate Merscope data (Fig. 4) are derived based on this R script. But the major pre-processing, inference, and post-processing steps shown in this R script are all consistent with those in our main spacia API.

```
Rscript [path/to/execute_spacia.R] \
	-i [path/to/input] \
	-r celltype1 \
	-s celltype2 \
	-g gene1 \
	-t [path/to/gene_cutoffs_prostate1.csv] \
	-o [path/to/output_celltype2-celltype1_gene1]
```
Use `-h` or `--help` to see detailed descriptions of options and inputs.

### Contact Us
If you have any suggestions/ideas for Spacia or are having issues trying to use it, please don't hesitate to reach out to us.

Yunguan Wang, yunguan[dot]wang@utsouthestern[dot]edu

Tao Wang, tao[dot]wang@utsouthestern[dot]edu
