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

### Quick start
Once the input data have been processed into the supported format, the full sprod workflow can be run by calling the `sprod.py` script. Sprod will first try to locate a single `.tif` image file in the input path. If one is found, sprod will assume it is the matching image and extract features from the image. If sprod cannot find an image, it will perform soft clustering on the spots and use the clustering probabilities as the input features for the denoising model. It evaluates interactions within the context of cell neighborhoods, where the ‘**receiver**’ cells are the cells of interest, and the cells from the neighborhood are referred to as "**sender**" cells. The **interactant** expressed in the receiver cells, through which the interactions are to be studied, are referred to as "**responder**", while the **interactant** expressed in the sender cells that potentially influence the responder genes are called signal “**sending**".

```
python [path/to/spacia_job.py] counts.txt cell_metadata.txt -rc celltype1 sc celltype2 -rf gene1 sf gene2
```

Here, `counts.txt` is a cell-by-gene matrix in TX format. We expect the data to be normliazed, if not, CPM normalization will be used.

`cell_metadata.txt` is a cell_by_metadata matrix in txt format in TXT format. Must contains `X` and `Y` columns for coordinates, and a `cell_type` columns, refering to the group designation of cells, is needed if '-rc' or '-sc' parameter are given.

`-rc` and `-sc` refer to **receiver** cells and **sender** cells, respectively.

`-rf` and `-sf` refer to **responder** and **sending** features. Here they are in forms of single genes. Spacia can also take pathways in the format of a list of genes as inputting features. 

### How Spacia defines interactant
Todo

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
