"""
Dependency
----------
R >= 4.0.2
    Rcpp, gcc
if running in biohpc, need to do module purge; module load shared first.
Python dependencies will be handeled automatically by installing the spacia package.

Test:
python [path/to/spacia.py] /project/shared/xiao_wang/projects/cell2cell_inter/data/Counts.txt \
    /project/shared/xiao_wang/projects/cell2cell_inter/data/Spot_metadata.txt \
    FGFR1 \
    -o, /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/spacia_py_test \
    -rc C_1 -sc C_2 \
    -sf FGF1,FGF2,FGF7,FGF9,FGF10,FGF11,FGF12,FGF13,FGF14,FGF17,FGF18 \
    # -cf /endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/spacia_py_test/input/input_cells.csv
"""
#%%
import os
import sys
import argparse
import logging
import csv
import json
from multiprocessing import Pool
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import scale
from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA

def spacia_worker(cmd):
    """
    worker function for multiprocessing.
    """
    os.system(cmd)
    # remove temp job input files
    # os.system('rm -f {}'.format(' '.join(spacia_job_inputs)))
    return

def cal_norm_dispersion(cts):
    '''
    Adapted from Scanpy _highly_variable_genes_single_batch.
    https://github.com/theislab/scanpy/blob/f7279f6342f1e4a340bae2a8d345c1c43b2097bb/scanpy/preprocessing/_highly_variable_genes.py
    '''
    mean, var = cts.mean(),  cts.std()
    mean[mean == 0] = 1e-12  # set entries equal to zero to small value
    dispersion = var / mean
    dispersion[dispersion == 0] = np.nan
    dispersion = np.log(dispersion)
    mean = np.log1p(mean)

    # all of the following quantities are "per-gene" here
    df = pd.DataFrame()
    df['means'] = mean
    df['dispersions'] = dispersion
    df['mean_bin'] = pd.cut(df['means'], bins=20)
    disp_grouped = df.groupby('mean_bin')['dispersions']
    disp_mean_bin = disp_grouped.mean()
    disp_std_bin = disp_grouped.std(ddof=1)
    # retrieve those genes that have nan std, these are the ones where
    # only a single gene fell in the bin and implicitly set them to have
    # a normalized disperion of 1
    one_gene_per_bin = disp_std_bin.isnull()
    disp_std_bin[one_gene_per_bin.values] = disp_mean_bin[
        one_gene_per_bin.values
    ].values
    disp_mean_bin[one_gene_per_bin.values] = 0
    # actually do the normalization
    df['dispersions_norm'] = (
        df['dispersions'].values  # use values here as index differs
        - disp_mean_bin[df['mean_bin'].values].values
    ) / disp_std_bin[df['mean_bin'].values].values
    return df.means.values, df.dispersions_norm.values

def calculate_neighbor_radius(
    spot_meta, sample_size=1000, target_n_neighbors=10, margin=10, r_min=1, r_max=1000
):
    r_samples = (
        spot_meta.loc[np.random.choice(spot_meta.index, sample_size, False)]
        .copy()
        .loc[:, ["X", "Y"]]
    )
    r_sample_dist = cdist(r_samples, spot_meta.loc[:, ["X", "Y"]])
    n_steps = 1
    while n_steps <= 100:
        r_next = (r_max + r_min) / 2
        n_neighbors = np.median(np.sum((r_sample_dist <= r_next), axis=1))
        if n_neighbors == target_n_neighbors:
            break
        elif (r_max - r_min) / 2 < margin:
            break
        n_steps += 1
        nn_r_min = np.median(np.sum((r_sample_dist <= r_min), axis=1))
        nn_r_max = np.median(np.sum((r_sample_dist <= r_max), axis=1))
        if (n_neighbors > target_n_neighbors) == (nn_r_max > target_n_neighbors):
            r_max = r_next
        elif (n_neighbors > target_n_neighbors) == (nn_r_min > target_n_neighbors):
            r_min = r_next
    return r_next


def find_sender_candidates(r_cells, s_cells, locations, dist_cutoff=30):
    """
    r_cells, s_cells: list of spot ids.
    locations: pd.DataFrame of X, Y locations
    """
    pip = pd.Series(dtype=object)
    n_chunks = int(np.ceil(len(r_cells)/10000))
    for c in range(n_chunks):
        r_cells_chunk = r_cells[c*10000 : (c+1)*10000]
        rl = locations.loc[r_cells_chunk].values.reshape(-1, 2)
        sl = locations.loc[s_cells].values.reshape(-1, 2)
        dist_to_senders = cdist(rl,sl)
        crit = dist_to_senders <= dist_cutoff
        r_candidates = r_cells_chunk[crit.any(axis=1)]
        crit = crit[crit.any(axis=1)]
        _pip = pd.Series([s_cells[crit[i]].tolist() for i in range(len(crit))])
        _pip.index = r_candidates
        pip = pd.concat([pip, _pip])
    return pip


def preprocessing_counts(
    counts, ntotal_cutoff=100, n_genes_cutoff=20, n_cells_cutoff=10
):
    """
    Simple QC based on total counts, num_genes expressed in cell and num of cells a gene is expressed.
    """
    counts = counts.T.groupby(counts.columns).mean().T
    # add the total counts per cell as observations-annotation to counts
    n_total = counts.sum(axis=1)
    n_genes = np.sum(counts > 0, axis=1)
    n_cells = np.sum(counts > 0, axis=0)
    keep_cells = np.zeros(shape=[counts.shape[0]], dtype=bool)
    keep_cells[:] = True
    keep_genes = np.zeros(shape=[counts.shape[1]], dtype=bool)
    keep_genes[:] = True
    for i, df, cutoff, qc_type in zip(
        [0, 0, 1],
        [n_total, n_genes, n_cells],
        [ntotal_cutoff, n_genes_cutoff, n_cells_cutoff],
        ["Total counts", "Total genes", "Number of cells expressed in"],
    ):
        crit = df < cutoff
        num_bad = crit.sum()
        if num_bad > 0:
            print(
                "{} {} are dropped because {} is less than {}.".format(
                    num_bad, "Cells" if i == 0 else "Genes", qc_type, cutoff
                )
            )
            if i == 0:
                keep_cells = np.logical_and(keep_cells, ~crit)
            else:
                keep_genes = ~crit.values
    filtered_counts = counts.loc[keep_cells, keep_genes]
    filtered_cpm = filtered_counts.apply(lambda x: 1e4 * x / x.sum(), axis=1)
    return filtered_cpm

def get_corr_agg_genes(corr_agg, cpm, cells, g, top_corr_genes, agg_method):
    if corr_agg:
        print('Constructing pathway using correlation aggregation')
        corr = 1-cdist(
            cpm.loc[cells, g].values.reshape(1, -1),
            cpm.loc[cells].T,
            metric="correlation",
        )[0]
        pathway_genes = pd.Series(corr, index = cpm.columns)
        if agg_method == 'simple':
            pathway_genes = pathway_genes[pathway_genes>0]
        else:
            pathway_genes = abs(pathway_genes)
        pathway_genes = pathway_genes.sort_values(
            ascending=False
            )[:top_corr_genes].index.tolist()
        pathway_name = g + "_correlated_genes"
    else:
        logging.warning("Correlation aggregation is turned off and this pathway has only one gene. This is not recommended.")
        pathway_name = g
        pathway_genes = [g]
    return pathway_genes, pathway_name

def contruct_pathways(
    cpm,
    receiver_candidates,
    sender_candidates,
    receiver_features,
    sender_features,
    agg_method,
    pca_gene = None
):
    receiver_pathways = {}
    sender_pathways = {}
    for pathway_dict, pathway_features, pathway_type in zip(
        [receiver_pathways, sender_pathways],
        [receiver_features, sender_features],
        ["Receiver", "Sender"],
    ):
        cells  = receiver_candidates if pathway_type == "Receiver" else sender_candidates
        if pathway_features == 'pca':
            pathway_exp = cpm.loc[cells,:]
            # Remove genes with all 0s
            pathway_exp = pathway_exp.T[pathway_exp.std() > 0].T
            # Calculate normalized dispersion and use it as cutoff
            pca = PCA(20).fit(scale(pathway_exp))
            pcc = pca.components_
            pcc = pd.DataFrame(
                pcc,
                index = ['PC_' + str(i+1) for i in range(pcc.shape[0])],
                columns = pathway_exp.columns
                )
            pcc = pcc.apply(lambda x: x/x.std())
            if pca_gene is not None:
                kept_pcs = abs(pcc.iloc[:5][pca_gene]).sort_values().index[:3].tolist()
            else:
                kept_pcs = pcc.index
            pcc = pcc.loc[kept_pcs]
            pathway_exp_pca = pd.DataFrame(
                np.matmul(scale(pathway_exp), pcc.T.values),
                index = pathway_exp.index,
                columns = pcc.index
            )
            pathway_dict[pathway_type + '_pc'] = pcc
            pathway_dict[pathway_type + '_y'] = pathway_exp_pca
            
        elif pathway_features is None:
            print(
                "{} features is not provided, use gene modules as pathways.".format(pathway_type)
            )
            # Get gene modules
            
            pathway_exp = cpm.loc[cells,:]
            # Remove genes with all 0s
            pathway_exp = pathway_exp.T[pathway_exp.std() > 0].T
            
            # Calculate normalized dispersion and use it as cutoff
            mean, ndisp = cal_norm_dispersion(pathway_exp)
            top_expressed_genes = (mean>=0.05) & (ndisp>0.05)

            pathway_exp = pathway_exp.loc[:,top_expressed_genes].copy()
            pathway_exp.loc[:,:] = scale(pathway_exp) # zscoring
            # correlation distance cutoff at 0.15
            gene_clusters = AgglomerativeClustering(
                None,
                affinity='correlation',
                linkage="complete", 
                distance_threshold=0.9,
                ).fit_predict(pathway_exp.T)
            
            # clean up clusters, removing singleton and big clusters
            vc = pd.Series(gene_clusters).value_counts()
            vc = vc.index[(vc>=5) & (vc<=100)]
            # assign genes not in a valid cluster to cluster -1
            gene_clusters = np.array([x if x in vc else -1 for x in gene_clusters])

            # construct sender_pathway
            n_c = np.unique(gene_clusters[gene_clusters!=-1]).shape[0]
            print(
                "Cosntruct {} {} pathways from gene modules".format(pathway_type,n_c)
                )
            for cluster in np.unique(gene_clusters):
                # ignore bad gene cluster -1
                if cluster == -1:
                    continue
                gene_mask = gene_clusters == cluster
                pathway_dict["module_" + str(cluster + 1)] = pathway_exp.columns[
                    gene_mask
                ].tolist()
        elif pathway_features[-4:] == ".csv":
            print("Cosntruct {} pathways from file...".format(pathway_type))
            with open(pathway_features) as csvfile:
                spamreader = csv.reader(csvfile, delimiter=",")
                for row in spamreader:
                    pathway_name = row[0]
                    pathway_genes = [x for x in row[1:] if x != ""]
                    pathway_genes = [x for x in pathway_genes if x in cpm.columns]
                    if len(pathway_genes) > 1:
                        pass
                    # If only one gene is present, will use correlations.
                    elif len(pathway_genes) == 1:
                        g = pathway_genes
                        if g not in cpm.columns:
                            print("{} not found in expression data.".format(g))
                            continue
                        pathway_genes, pathway_name = get_corr_agg_genes(
                            corr_agg, cpm, cells, g, top_corr_genes, agg_method)
                    else:
                        continue # just to handle blank lines
                    pathway_dict[pathway_name] = pathway_genes
        elif "|" in pathway_features:
            print("Cosntruct 1 {} pathway from input genes".format(pathway_type))
            genes = pathway_features.split("|")
            pathway_dict[pathway_type + "_pathway"] = genes
        else:
            print("Cosntruct {} pathways from each input gene".format(pathway_type))
            genes = pathway_features.split(",")
            for g in genes:
                if g not in cpm.columns:
                    print("{} not found in expression data.".format(g))
                    continue
                pathway_genes, pathway_name = get_corr_agg_genes(
                    corr_agg, cpm, cells, g, top_corr_genes, agg_method)
                pathway_dict[pathway_name] = pathway_genes
    return receiver_pathways, sender_pathways


class StreamToLogger(object):
    """
    Fake file-like stream object that redirects writes to a logger instance.
    """

    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ""

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())

    def flush(self):
        pass

def remove_outliers(betas):
    outlier_rows = []
    for col in betas:
        cutoff_l = betas[col].mean() - 5* betas[col].std()
        cutoff_2 = betas[col].mean() + 5* betas[col].std()
        outlier_rows +=betas.index[
            (betas[col]>cutoff_2) | (betas[col]<cutoff_l)
            ].tolist()
    outlier_rows = list(set(outlier_rows))
    betas = betas[~betas.index.isin(outlier_rows)]
    return betas
#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Main function for running spacia. There are several running modes and spacia will \
            automatically decide which one to use based on the input parameters. Below are the four modes: \n\t\
            (1) No aggregation: The user either provides one or several single genes, \
            for sending/receiving pathways. These will be used in the analysis directly; \
            These genes have to be positively correlated. \n\t\
            (2) Correlation-driven aggregation: The users will provide a list of single genes. \
            Each single gene will be expanded, in spacia internally, to a pathway, \
            based on highly positively correlated genes. The user will provide a correlation cutoff \
            to select the highly positively correlated genes for each seed gene; \n\t\
            (3) Knowledge-driven aggregation: The users can provide a .csv file, \
            detailing what genes to group into pathways, for sending/receiving pathways. \
            Users can refer to CellChat, CellphoneDB and other appropriate resources to collect \
            such pre-defined pathways and their member genes. But users need to take care to only \
            include genes that are supposedly positively correlated in the same pathway; \n\t\
            (4) Clustering-driven aggregation: If the users do not provide any information, \n\t\
            spacia will perform un-supervised clustering and divide all genes into modules based on \
            positive correlation, based on the positive correlation cutoff the users provided. \
            Each module will be treated as a pathway.",
    )

    parser.add_argument(
        "counts",
        help="Path for gene expression data, spots by genes, must be normlized. \
            TXT format",
    )

    parser.add_argument(
        "spot_meta",
        help="Path for spot positional information, spots by features. TXT format. \
            must have 'X', 'Y' columns for coordinates. 'cell_type' columns is needed \
            if running with '-rc' and '-sc' parameters. 'cell_type' refers to the group designation of cells, \
            e.g., type of cells. The user can specify which group of cells to use for spacia"
    )

    parser.add_argument(
        "--receiver_features",
        "-rf",
        type=str,
        default=None,
        help="Receiver gene(s). Can be: \
            1) a single gene. If 'corr_agg' is turned off, spacia will run in mode (1); \
            otherwise mode (2); \
            2) multiple genes, sep by ',' each is a seed of one receiver pathway. \
            Spacia will run in mode (2); \
            3) multiple genes, sep by '|' used together as one pathway. Spacia will \
            run in mode (3); \
            4) a csv file each row for a receiver pathway and columns are genes, \
            first column contains pathways names. Spacia will run in mode (3) if \
            the pathway contains multiple genes. If there is only one gene in the pathway,\
            spacia will run in mode (1) if 'corr_agg' is turned off otherwise in mode (2); \
            4) 'all', spacia will use all genes in the count matrix; \
            6) if nothing is passed, spacia will run in mode (4).",
    )

    parser.add_argument(
        "--sender_features",
        "-sf",
        type=str,
        default=None,
        help="Sender gene(s). At least two pathways are recommended. Can be: \
            1) a single gene. If 'corr_agg' is turned off, spacia will run in mode (1); \
            otherwise mode (2); \
            2) multiple genes, sep by ',' each is a seed of one receiver pathway. \
            Spacia will run in mode (2); \
            3) multiple genes, sep by '|' used together as one pathway. Spacia will \
            run in mode (3); \
            4) a csv file each row for a receiver pathway and columns are genes, \
            the first column contains pathway names. Spacia will run in mode (3) if \
            the pathway contains multiple genes. If there is only one gene in the pathway,\
            spacia will run in mode (1) if 'corr_agg' is turned off otherwise in mode (2); \
            5) if 'pca' is provided as input, use transform data using top 20 principle \
            components and use each component as a sender pathway. \
            6) if nothing is passed, spacia will run in mode (4).",
    )

    parser.add_argument (
        "--pca_gene",
        "-pg",
        help = "Additional gene to be included in the pca mode. Note this will cause spacia \
            to use only the 3 top pcs."
    )
    parser.add_argument(
        "--output_path", "-o", type=str, default="spacia", help="Output path"
    )

    parser.add_argument(
        "--dist_cutoff",
        "-d",
        type=float,
        default=None,
        help="Distance cutoff for finding potential interacting cells.",
    )

    parser.add_argument(
        "--receiver_exp_cutoff",
        "-rec",
        type=str,
        default=0.5,
        help="Quantile cutoff used to threshold receiver expression. If 'auto', a bimodel \
            distribution will be fitted to find the cutoff.",
    )

    parser.add_argument(
        "--num_corr_genes",
        "-nc",
        type=int,
        default=100,
        help="Number of correlated gene to use in calculating receiver pathway expression.\
            This option only matters if the user has not toggle off the 'corr_agg' option.",
    )

    parser.add_argument(
        "--n_neighbors",
        "-n",
        type=float,
        default=10,
        help="Expected number of nearest neighbors for most of the cells. The \
        exact number could vary from cell to cell.",
    )

    parser.add_argument(
        "--mcmc_params",
        "-m",
        type=str,
        default="50000,25000,10,3",
        help="MCMC parameters, four values packed here are {ntotal,nwarm,nthin,nchain}",
    )


    parser.add_argument(
        "--receiver_cluster",
        "-rc",
        type=str,
        default=None,
        help="Name of receiver cell_type, must be in spot_metadata.",
    )

    parser.add_argument(
        "--sender_cluster",
        "-sc",
        type=str,
        default=None,
        help="Name of sender cell_type, must be in spot_metadata.",
    )
    
    parser.add_argument(
        "--bag_size",
        "-b",
        type=int,
        default=2,
        help="Minimal bag size for sender cells.",
    )
    
    parser.add_argument(
        "--number_bags",
        "-nb",
        type=int,
        default=5000,
        help="Number of bags used in multiple instance learning.",
    )

    parser.add_argument(
        "--cellid_file",
        "-cf",
        type=str,
        default=None,
        help="Name of a csv file for receiver and sender cell ids. \
             first columns for receiver cells, second for sender cells. Will be \
             overridden if 'receiver_cluster' and 'sender_cluster' are given",
    )

    parser.add_argument(
        "--keep_intermediate",
        "-k",
        action="store_false",
        default=True,
        help="Whether the model_input folder should be deleted.",
    )

    parser.add_argument(
        "--corr_agg",
        "-ca",
        action="store_false",
        default=True,
        help="This is a toggle to turn off correlation based aggregation. If it is \
            turned off, spacia evaluate the effect of one gene in the sender cells on \
            another gene in the receiver cells.",
    )
    
    parser.add_argument(
        "--corr_agg_method",
        "-cm",
        default='simple', # simple or weighted
        help="How the receiver gene expression will be aggregated, if 'simple', top positively \
            correlated genes will be averages. If 'weighted', top correlated genes by absolute \
            correlation will be summed into the seed genes expression using weighted averages",
    )

    parser.add_argument (
        "--plot_mcmc",
        action = "store_true",
        default = False,
        help = "Optional argument for plotting b and beta's trace plots, density plots, \
         autocorrelation plots, and PSRF plots."
    )
    
    parser.add_argument (
        "--debug_plots",
        action = "store_true",
        default = False,
        help = "Optional argument for making debug plots for exp_receiver cutoff selection."
    )

    parser.add_argument (
        "--ext",
        type = str,
        default = 'pdf',
        help = "File formats for the mcmc plots to be saved.Can either be a device function \
         (e.g. png), or one of eps, ps, tex (pictex), pdf, jpeg, tiff, png, bmp, svg or wmf (windows only)"
    )

    ######## Setting up ########
    # Debug param
    
#%%
    ######## Setting up ########
    args = parser.parse_args()
    # args = parser.parse_args(
    #     [
    #         '/home2/s190548/work_personal/software/cell2cell_inter/code/test/input/Counts.txt',
    #         '/home2/s190548/work_personal/software/cell2cell_inter/code/test/input/spacia_metadata.txt',
    #         '-rc', 'A',
    #         '-sc', 'B',
    #         '-rf', 'gene1',
    #         '-sf', 'gene2,gene3',
    #         '-d', '5',
    #         '-m', '5000,2500,10,2',
    #         '-nc', '20',
    #         '-o', '/home2/s190548/work_personal/software/cell2cell_inter/code/test/test_output'
    #         ])
    counts = args.counts
    spot_meta = args.spot_meta
    receiver_cluster = args.receiver_cluster
    sender_cluster = args.sender_cluster
    cellid_file = args.cellid_file
    # pathway_lib = args.pathway_lib
    output_path = args.output_path
    n_neighbors = args.n_neighbors
    receiver_features = args.receiver_features
    sender_features = args.sender_features
    top_corr_genes = args.num_corr_genes
    dist_cutoff = args.dist_cutoff
    receiver_exp_cutoff = args.receiver_exp_cutoff
    receiver_exp_cutoff = receiver_exp_cutoff if receiver_exp_cutoff == 'auto' else float(receiver_exp_cutoff)
    mcmc_params = args.mcmc_params
    corr_agg = args.corr_agg
    ntotal, nwarm, nthin, nchain = mcmc_params.split(",")
    keep = args.keep_intermediate
    plot_mcmc = 'T' if args.plot_mcmc else 'F'
    ext = args.ext
    plot_debug = args.debug_plots
    corr_agg_method = args.corr_agg_method
    bag_size = args.bag_size
    nb = args.number_bags
    pca_gene = args.pca_gene
    np.random.seed(0)

    # Checking inputs
    assert corr_agg_method in ['simple','weighted'], \
        f"'corr_agg_method' must be either 'simple' or 'weighted'!"  
    
    intermediate_folder = os.path.join(output_path, "model_input")
    if not os.path.exists(intermediate_folder):
        os.makedirs(intermediate_folder)

    # Setting up logs
    log_fn = os.path.join(output_path, "spacia_log.txt")
    if os.path.exists(log_fn):
        os.remove(log_fn)
    logging.basicConfig(
        filename=log_fn,
        format="%(asctime)s,%(levelname)s:::%(message)s",
        datefmt="%H:%M:%S",
        level="INFO",
    )
    # print(args)

    # getting script path for supporting codes.
    spacia_path = os.path.abspath(__file__)
    wrapup_script = spacia_path.replace("spacia.py", "wrapup_spacia_results.py")
    spacia_path = "/".join(spacia_path.split("/")[:-1]) + "/spacia"
    spacia_script = os.path.join(spacia_path, "spacia_job.R")

    # redirects stdout and stderr to logger
    stdout_logger = logging.getLogger("STDOUT")
    sl = StreamToLogger(stdout_logger, logging.INFO)
    sys.stdout = sl
    stderr_logger = logging.getLogger("STDERR")
    sl = StreamToLogger(stderr_logger, logging.ERROR)
    sys.stderr = sl

    ######## Processing counts and receiver and sender cells ########
    # Processing counts and spot_metadata
    print('Processing expression counts.')
    counts = pd.read_csv(counts, index_col=0, sep="\t")
    spot_meta = pd.read_csv(spot_meta, index_col=0, sep="\t")
    if not all(x in spot_meta.columns for x in ['X','Y','cell_type']):
        raise ValueError(
            "Metadata must have ['X','Y','cell_type'] columns!"
        )
    # TODO: added a tag to allow normalization
    if counts.max().max() > 1000:
        # cpm = preprocessing_counts(counts)
        UserWarning(
            'input gene expression data does not seem in log1cpm format'
            )
    else:
        cpm = counts
    cpm, spot_meta = cpm.align(spot_meta, join="inner", axis=0)

    # find candidate receiver and sender cells
    if dist_cutoff is None:
        dist_cutoff = calculate_neighbor_radius(
            spot_meta.iloc[:, :2], target_n_neighbors=n_neighbors
        )
        print(
            "Maximal distance for {} expected neighbors is {:.2f}".format(
                n_neighbors, dist_cutoff
            )
    )
    # catch error where a wrong cell cluster name is provided.
    for c_name in [receiver_cluster, sender_cluster]:
        if c_name not in spot_meta.cell_type.unique():
            raise ValueError('{} not found in cell types!'.format(spot_meta))
        
    if (sender_cluster is not None) & (sender_cluster is not None):
        r_cells = spot_meta[spot_meta.cell_type == receiver_cluster].index
        s_cells = spot_meta[spot_meta.cell_type == sender_cluster].index
    elif cellid_file is not None:
        cellids = pd.read_csv(cellid_file, header=None)
        r_cells = cellids.iloc[:, 0].dropna().values
        s_cells = cellids.iloc[:, 1].dropna().values
    else:
        raise ValueError(
            "Must provide both receiver and sender clusters, or a file with their ids."
        )

    r2s_matrix = find_sender_candidates(
        r_cells, s_cells, spot_meta[["X", "Y"]], dist_cutoff
    )
    receiver_cell_for_cutoff = r2s_matrix.index.tolist()
    print('Limiting bags to those with at least {} sender cells'.format(bag_size))
    r2s_matrix = r2s_matrix[r2s_matrix.apply(len) >= bag_size]
    if r2s_matrix.shape[0] < 500:
        raise ValueError('Number of total bags is too small, job killed.')
    elif r2s_matrix.shape[0]> nb:
        print('Subsample bags for Spacia.')
        r2s_matrix = r2s_matrix.sample(nb, replace=False)
    sender_candidates = list(set(r2s_matrix.sum()))
    receiver_candidates = r2s_matrix.index.tolist()

    ######## Preparing spacia_job.R inputs ########
    # Contruct sender and receiver pathways
    if receiver_features == 'all':
        receiver_features = ','.join(cpm.columns)
    receiver_pathways, sender_pathways = contruct_pathways(
        cpm, 
        receiver_candidates, 
        sender_candidates, 
        receiver_features, 
        sender_features,
        corr_agg_method,
        pca_gene
    )
    # If no receiver pathways are found, abort.
    if len(receiver_pathways.keys()) == 0:
        print('None of the genes in the provided receiver pathways are found in \
            the expression matrix, please modify the input and try again.')
        raise ValueError()
        
    print('Writing spacia_job.R inputs to the model_input folder.')
    dist_sender_fn = os.path.join(intermediate_folder, "dist_sender.json")
    metadata_fn = os.path.join(intermediate_folder, "metadata.txt")
    exp_sender_fn = os.path.join(intermediate_folder, "exp_sender.json")
    # Calculate each receiver sender pair distances
    dist_r2s = r2s_matrix.to_frame().apply(
        lambda x: (cdist(
            spot_meta.loc[x.name, :"Y"].values.reshape(-1, 2),
            spot_meta.loc[x[0], :"Y"].values.reshape(-1, 2),
        )[0]/dist_cutoff).round(5), # normalize distance to 0-1
        axis=1,
    )
    sender_dist_dict = {}
    for i in dist_r2s.index:
        sender_dist_dict[i] = dist_r2s[i].tolist()

    # contruct and save metadata
    meta_data = spot_meta.loc[receiver_candidates, :"Y"]
    meta_data["Sender_cells"] = r2s_matrix.loc[receiver_candidates].apply(",".join)
    meta_data_senders = spot_meta.loc[sender_candidates, :"Y"]
    meta_data = meta_data.append(meta_data_senders)

    # contruct and save sender exp
    if sender_features == 'pca':
        sender_pathway_exp = sender_pathways['Sender_y']
        if pca_gene is not None:
            sender_pathway_exp[pca_gene] = cpm.loc[sender_pathway_exp.index, pca_gene]
        sender_pathway_exp.loc[:,:] = scale(sender_pathway_exp)
    else:
        sender_pathway_exp = pd.DataFrame(
            index=sender_candidates, columns=sender_pathways.keys()
        )
        for key in sender_pathway_exp.columns:
            sender_pathway_exp[key] = scale(
                cpm.loc[sender_candidates, sender_pathways[key]].mean(axis=1)
            )
        
    # # Add one dummy pathway as control
    # dummy_pathway = np.random.normal(
    #     scale=0.01,
    #     size=sender_pathway_exp.shape[0]
    # )
    # sender_pathway_exp['dummy'] = dummy_pathway
        
    sender_exp = (
        r2s_matrix.to_frame()
        .apply(lambda x: sender_pathway_exp.loc[x[0],].values.round(5).tolist(), axis=1)
        .to_dict()
    )
    
    ######## Write spacia_job.R jobs ########
    # construct receiver expression and the job commands
    spacia_jobs = []
    spacia_job_folders = []
    for rp in receiver_pathways.keys():
        job_id = rp
        job_folder = os.path.join(output_path, job_id)
        spacia_job_folders.append(job_folder)
        
        # Check if the current rp is already done
        log_path = os.path.join(job_folder, job_id + '_log.txt')
        job_finished = False
        if os.path.exists(log_path):
            with open(log_path, 'r') as f:
                log = f.readlines()
                job_finished = any(
                    list(map(lambda x: 'Time difference' in x, log)))
        if job_finished:
            print(job_id + ' is already finished and will be skipped.')
            continue
        
        exp_receiver_fn = os.path.join(
            intermediate_folder, job_id + "_exp_receiver.csv"
        )
        # Getting receiver exp
        rp_genes = receiver_pathways[rp]
        # aggregate gene expression
        if corr_agg_method == 'simple':
            receiver_exp = cpm.loc[receiver_cell_for_cutoff, rp_genes].mean(axis=1)
        else:
            corr = cpm.loc[receiver_cell_for_cutoff, rp_genes].corr()[rp.split('_')[0]]
            receiver_exp = np.matmul(
                cpm.loc[receiver_cell_for_cutoff, rp_genes],corr
                )
        # Decide receiver exp cutoff
        # Debug codes
        # print(receiver_exp.head())
        # print(receiver_exp_cutoff)
        rf_to_drop = []
        if receiver_exp_cutoff == 'auto':
            print(
                'Estimating {} expression cutoff by fitting a bimodal distribution...'.format(rp)
                )
            gm = GaussianMixture(n_components=2, random_state=0).fit(
                receiver_exp.values.reshape(-1,1))
            labels = gm.predict(receiver_exp.values.reshape(-1,1))
            # check bimodality and calculate cutoff
            sd1 = receiver_exp[labels==0].std()
            sd2 = receiver_exp[labels==1].std()
            m1, m2 = gm.means_.flatten()
            if m1 > m2:
                m1, m2 = m2, m1
                sd1, sd2 = sd2, sd1
            if m2-1*sd2 <= m1+1*sd1:
                # If not bimodal, use median
                print('{} expression is likely not bimodal!'.format(rp))
                print('Using m1 + 1sd cutoff value.')
                cutoff = m1+1*sd1
            else:
                cutoff = (m1+m2)/2
                
            # For pathways whose expression are very expreme, use median as cutoff
            if (
                (labels.sum() > 0.9 * receiver_exp.shape[0]) or 
                (labels.sum() < 0.1 * receiver_exp.shape[0])
            ):
                # print('Receiver expression too extreme, job skipped')
                print('Receiver expression maybe too extreme.')
                # rf_to_drop.append(rp)
                cutoff = receiver_exp.quantile(0.5)
        else:
            cutoff = receiver_exp.quantile(receiver_exp_cutoff)
    
        if plot_debug:
            receiver_exp.hist(bins=20,density=True)
            plt.plot((cutoff,cutoff), (0,2))
            plt.savefig(
                os.path.join(intermediate_folder, job_id + "_exp_receiver_dist.pdf"))
            plt.close()
            
        receiver_exp = receiver_exp > cutoff
        receiver_exp = receiver_exp + 0
        receiver_exp = receiver_exp[receiver_candidates]
        receiver_exp.to_csv(exp_receiver_fn, header=None, index=None)

        spacia_output_path = os.path.join(output_path, job_id)
        if not os.path.exists(spacia_output_path):
            os.makedirs(spacia_output_path)
        spacia_jobs.append(
            " ".join(
                [
                    "Rscript",
                    spacia_script,
                    spacia_path + "/",
                    exp_sender_fn,
                    dist_sender_fn,
                    exp_receiver_fn,
                    job_id,
                    str(ntotal),
                    str(nwarm),
                    str(nthin),
                    str(nchain),
                    spacia_output_path + "/",
                    plot_mcmc,
                    ext,
                ]
            )
        )
    
    with open(os.path.join(output_path, 'spacia_r.log'), 'w') as f:
        f.write('\n'.join(spacia_jobs)) # Save the actual jobs for debug purpose
        
    # Save receiver and sender pathways for reference
    # remove receiver genes from receiver pathway
    # for key in rf_to_drop:
    #     del receiver_pathways[key]
    # sender_pathways['dummy'] = [] # add dummy pathway
    for pathway_dict, fn in zip(
        [receiver_pathways, sender_pathways],
        ["receiver_pathways.json", "sender_pathways.json"],
    ):
        if (fn == "sender_pathways.json") & (sender_features == 'pca'):
            pc_fn = os.path.join(intermediate_folder, 'sender_pc.csv')
            sender_pathways['Sender_pc'].to_csv(pc_fn)
            # prepare dummy sender_pathway.json for pca mode
            pathway_dict = pd.DataFrame(
                index=sender_pathways['Sender_pc'].index.tolist())
            pathway_dict = pathway_dict.to_dict(orient='index')
            if pca_gene is not None:
                pathway_dict[pca_gene] = {}

        with open(os.path.join(intermediate_folder, fn), "w") as fp:
            json.dump(pathway_dict, fp)
            
    # Writing spacia R job inputs common for each receiver pathways
    # job metadata
    meta_data.to_csv(metadata_fn, sep='\t')
    
    # sender distance and expression json (list of lists)
    with open(dist_sender_fn, "w") as fp:
        json.dump(sender_dist_dict, fp)
        
    with open(exp_sender_fn, "w") as fp:
        json.dump(sender_exp, fp)
    
    ######## Proceed with spacia_job.R ########
    # Run all spacia R jobs
    print('Running spacia_R MCMC MIL models.')
    with Pool(16) as p:
        _ = p.map(spacia_worker, spacia_jobs)
    
    ######## Collect all results ########
    print('Collecting results.')
    meta_data = pd.read_csv(metadata_fn, index_col=0, sep="\t")
    with open(os.path.join(intermediate_folder, "sender_pathways.json"), "r") as f:
        sender_pathways_names = json.load(f).keys()

    interactions_template = (
        meta_data.dropna(subset=["Sender_cells"])
        .Sender_cells.str.split(",", expand=True)
        .stack()
        .reset_index()
    )
    interactions_template.columns = ["Receiver", "x", "Sender"]
    interactions_template = interactions_template[["Receiver", "Sender"]]

    print('Spacia_R_results at: \n\t{}'.format('\n\t'.join(spacia_job_folders)))
    pathways = pd.DataFrame()
    interactions = pd.DataFrame()
    b_plus_fdr = pd.DataFrame()
    for fd in spacia_job_folders:
        job_id = fd.split('/')[-1]
        # aggregating beta for different receiver pathways
        try:
            res_beta = pd.read_csv(os.path.join(fd, job_id + "_beta.txt"), sep="\t")
            res_beta = remove_outliers(res_beta)
            res_beta = res_beta.apply(lambda x: x/x.std()).mean()
        except:
            print('{} failed without outputs!'.format(job_id))
            continue
        res_beta = res_beta.reset_index()
        res_beta.index = [job_id] * res_beta.shape[0]
        res_beta.columns = ["Sender_pathway", "Beta"]
        
        # no longer needed
        # # Fix an issue where there is only one sender pathway and the sender_exp
        # # has a dummy pathway expression.
        # if len(sender_pathways_names) == 1:
        #     sender_pathways_names = list(sender_pathways_names) + ['dummy']
            
        res_beta.Sender_pathway = sender_pathways_names
        pathways = pathways.append(res_beta)

        # aggregating primamy instances for different receiver pathways
        pip_res = pd.read_csv(
            os.path.join(fd, job_id + "_pip.txt"), sep="\t"
        ).mean(axis=1)
        assert (
            pip_res.shape[0] == interactions_template.shape[0]
        ), "Spaca results don't match input!"
        _interactions = interactions_template.copy()
        _interactions.index = [job_id] * _interactions.shape[0]
        _interactions["Primary_instance_score"] = pip_res.values
        interactions = interactions.append(_interactions)

        # aggregating b and FDR for different receiver pathways
        pred_b = (
            pd.read_csv(os.path.join(fd, job_id + "_b.txt"), sep="\t")
            .iloc[:, 1]
            .mean()
        )
        fdr = pd.read_csv(os.path.join(fd, job_id + "_FDRs.txt"), sep="\t")
        fdr = fdr.reset_index()
        fdr.index = [job_id] * fdr.shape[0]
        fdr.columns = ["Theta_cutoff", "FDR"]
        fdr.Theta_cutoff = fdr.Theta_cutoff / 10
        fdr["b"] = pred_b
        b_plus_fdr = b_plus_fdr.append(fdr)

        pathways.to_csv(os.path.join(output_path, "Pathway_betas.csv"))
        interactions.to_csv(os.path.join(output_path, "Interactions.csv"))
        b_plus_fdr.to_csv(os.path.join(output_path, "B_and_FDR.csv"))
    
    # Remove model_input files
    if not keep:
        os.system("rm -rf {}".format(intermediate_folder))
# %%
# import matplotlib.pyplot as plt
# test = cpm.loc[receiver_candidates]
# genes = test.mean().sort_values(ascending=False).index[:50]
# _, ax = plt.subplots(10,5, figsize = (20,40))
# ax = ax.ravel()
# for i, g in enumerate(genes):
#     ax[i].hist(test[g], bins=30)
#     ax[i].set_title(g)
# # %%
