#%%
import os
import glob
import hashlib
import subprocess
import shlex
import sys
import argparse
import logging
import csv
import json
from concurrent.futures import ThreadPoolExecutor
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import scale
from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA
from scipy import stats
import pprint

def spacia_worker(cmd):
    """
    worker function for multiprocessing.
    """
    subprocess.run(shlex.split(cmd) if isinstance(cmd, str) else cmd, check=True)
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
    spot_meta, r_cells, s_cells, sample_size=1000, target_n_neighbors=10,
):
    r_spot_meta = spot_meta.loc[r_cells]
    s_spot_meta = spot_meta.loc[s_cells]
    if sample_size >= r_spot_meta.shape[0]:
        r_samples = spot_meta
    else:
        r_samples = (
            r_spot_meta.loc[np.random.choice(r_cells, sample_size, False)]
            .copy()
            .loc[:, ["X", "Y"]]
        )
    r_sample_dist = cdist(r_samples, s_spot_meta.loc[:, ["X", "Y"]])
    r_min = np.min(r_sample_dist)
    r_max = np.max(r_sample_dist)
    margin = np.std(r_sample_dist)*0.01
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
    # Radius queries avoid allocating receiver-by-sender distance matrices.
    # Sort sender indices to preserve the original column order and RNG mapping.
    sender_xy = locations.loc[s_cells, ["X", "Y"]].to_numpy()
    receiver_xy = locations.loc[r_cells, ["X", "Y"]].to_numpy()
    tree = cKDTree(sender_xy)
    neighbors = tree.query_ball_point(receiver_xy, dist_cutoff, return_sorted=True)
    sender_ids = np.asarray(s_cells)
    return pd.Series(
        [sender_ids[ids].tolist() for ids in neighbors if len(ids)],
        index=[cell for cell, ids in zip(r_cells, neighbors) if len(ids)],
        dtype=object,
    )


def load_detection_counts(path, sender_gene, receiver_gene):
    """Read only the two detection genes from a raw cell-by-gene TSV."""
    with open(path, newline='') as handle:
        header = next(csv.reader(handle, delimiter='\t'))
    genes = list(dict.fromkeys([sender_gene, receiver_gene]))
    if any(header.count(gene) != 1 or header.index(gene) == 0 for gene in genes):
        raise ValueError("Detection counts must contain each selected gene exactly once as a gene column")
    columns = sorted({0, *(header.index(gene) for gene in genes)})
    raw = pd.read_csv(path, sep='\t', index_col=0, usecols=columns)
    if not raw.index.is_unique or raw.index.hasnans:
        raise ValueError("Detection-count cell IDs must be unique and non-missing")
    return raw


def filter_sender_pairs(r2s_matrix, zero_policy='keep', detection_counts=None,
                        sender_gene=None, receiver_gene=None, min_senders=2):
    """Filter raw-count double-zero pairs, without modifying feature fitting.

    Returns (retained neighborhoods, pair audit, receiver audit). A pair removed
    because its remaining neighborhood is too small is audited separately.
    """
    if zero_policy not in ('keep', 'exclude_double_zero'):
        raise ValueError("Unknown zero policy: " + str(zero_policy))
    if zero_policy == 'keep':
        return r2s_matrix.copy(), pd.DataFrame(), pd.DataFrame()
    if detection_counts is None or not sender_gene or not receiver_gene:
        raise ValueError("Double-zero exclusion requires raw counts and one sender/receiver gene")
    if min_senders < 1:
        raise ValueError("Minimum sender count must be positive")
    if not detection_counts.index.is_unique or detection_counts.index.hasnans:
        raise ValueError("Detection-count cell IDs must be unique and non-missing")
    if not detection_counts.columns.is_unique:
        raise ValueError("Detection-count gene columns must be unique")
    genes = list(dict.fromkeys([sender_gene, receiver_gene]))
    if any(g not in detection_counts.columns for g in genes):
        raise ValueError("Selected detection gene is missing from raw counts")
    cells = pd.Index(list(dict.fromkeys(
        list(r2s_matrix.index) + [s for senders in r2s_matrix for s in senders])))
    missing = cells.difference(detection_counts.index)
    if len(missing):
        raise ValueError("Raw detection counts are missing cells: " + ', '.join(map(str, missing[:10])))
    raw = detection_counts.loc[cells, genes].apply(pd.to_numeric, errors='raise')
    values = raw.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any():
        raise ValueError("Detection counts must be finite, nonnegative integer transcript counts; missing is not zero")
    retained = {}
    pair_rows, receiver_rows = [], []
    for receiver, senders in r2s_matrix.items():
        receiver_count = raw.at[receiver, receiver_gene]
        sender_counts = raw.loc[senders, sender_gene].to_numpy()
        double_zero = (sender_counts == 0) & (receiver_count == 0)
        remaining = [s for s, drop in zip(senders, double_zero) if not drop]
        keep_receiver = len(remaining) >= min_senders
        if keep_receiver:
            retained[receiver] = remaining
        receiver_rows.append(dict(receiver=receiver, receiver_count=receiver_count,
                                  senders_before=len(senders), double_zero_pairs=int(double_zero.sum()),
                                  senders_after_pair_filter=len(remaining), retained=keep_receiver))
        for sender, count, drop in zip(senders, sender_counts, double_zero):
            action = ('excluded_double_zero' if drop else
                      'retained' if keep_receiver else 'excluded_small_neighborhood')
            pair_rows.append(dict(receiver=receiver, sender=sender, receiver_count=receiver_count,
                                  sender_count=count, double_zero=bool(drop), action=action))
    return (pd.Series(retained, dtype=object), pd.DataFrame(pair_rows),
            pd.DataFrame(receiver_rows))


def prepare_zero_policy(args, parser):
    """Validate explicit raw-count semantics and prevent incompatible resume."""
    config = dict(mode=args.zero_policy)
    if args.zero_policy == 'exclude_double_zero':
        for name, gene in [('sender', args.sender_features), ('receiver', args.receiver_features)]:
            if (not gene or gene in ('all', 'pca') or ',' in gene or '|' in gene
                    or gene.lower().endswith('.csv')):
                parser.error("exclude_double_zero requires one gene in -sf and one gene in -rf")
        if not args.detection_counts:
            parser.error("exclude_double_zero requires --detection-counts RAW_COUNTS.tsv")
        if args.bag_size < 1:
            parser.error("--bag_size must be positive")
        raw = load_detection_counts(args.detection_counts, args.sender_features, args.receiver_features)
        digest = hashlib.sha256()
        with open(args.detection_counts, 'rb') as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b''):
                digest.update(chunk)
        config.update(sender_gene=args.sender_features, receiver_gene=args.receiver_features,
                      detection_counts=os.path.realpath(args.detection_counts),
                      detection_counts_sha256=digest.hexdigest(), min_senders=args.bag_size)
    else:
        raw = None
    marker = os.path.join(args.output_path, 'zero_policy.json')
    if os.path.exists(marker):
        with open(marker) as handle:
            if json.load(handle) != config:
                raise ValueError("Zero policy, detection genes or raw counts differ from this output directory; use a fresh --output_path")
    elif args.zero_policy != 'keep' and glob.glob(os.path.join(args.output_path, '*', '*_log.txt')):
        raise ValueError("Existing results have no zero-policy record; use a fresh --output_path for exclusion")
    os.makedirs(args.output_path, exist_ok=True)
    with open(marker, 'w') as handle:
        json.dump(config, handle, indent=2)
    return raw, config


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
    n_pc = 20,
    pca_gene = None,
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
            pca = PCA(n_pc).fit(scale(pathway_exp))
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

def format_json(dict):
    f = pprint.pformat(dict, sort_dicts=False)
    f = f.replace('\'','"')
    return f

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

def process_b(df_b, spacia_res_path, chain_size, n_chains):
    df_b = df_b.groupby(df_b.index).first()
    remove = [x*chain_size for x in range(n_chains)]
    indiv_results = df_b.index.unique()
    planned = os.listdir(spacia_res_path)
    for fn in [
        'Interactions.csv', 'B_and_FDR.csv', 'spacia_log.txt', 
        'Pathway_betas.csv', 'spacia_r.log', 'model_input', 'zero_policy.json',
        'zero_policy_pairs.csv', 'zero_policy_receivers.csv', 'zero_policy_summary.json']:
        try:
            planned.remove(fn)
        except:
            continue
    if len([x for x in planned if x not in indiv_results]) > 0:
        print('Warning!! ',
            [x for x in planned if x not in indiv_results],
            ' are not found in results!!!')
    for fn in indiv_results:
        fn_b = os.path.join(spacia_res_path, fn, fn + '_b.txt')
        if not os.path.exists(fn_b):
            print(fn_b, ' is not found!!')
            continue
        # print(fn_b)
        indiv_b = pd.read_csv(fn_b, header=None, sep='\t', skiprows=1).iloc[:,2]
        indiv_b = indiv_b[~indiv_b.index.isin(remove)]
        arr = []
        for i in range(n_chains):
            chain_df = indiv_b[chain_size*i:chain_size*(i+1)]
            arr += chain_df.sample(50, replace=False).tolist()
        arr = np.array(arr)
        pval = stats.ttest_1samp(arr, 0, alternative='less')[1]
        df_b.loc[fn, 'pval'] = pval
    df_b['pval_adj'] = p_adjust_bh(df_b['pval'])
    return df_b

def process_beta(
    pathway_beta, spacia_res_path, chain_size, n_chains, mode = 'pca'):
    indiv_results = pathway_beta.index.unique()
    remove = [x*chain_size for x in range(n_chains)]
    if mode == 'pca':
        pathway_beta = {
            'RG':[], 
            'Sender_pathway': [], 
            'Beta': [],
            'pval': [],
        }
    for rg in indiv_results:
        fn_beta = os.path.join(spacia_res_path, rg, rg + '_beta.txt')
        if not os.path.exists(fn_beta):
            continue
        df_beta = pd.read_csv(fn_beta, sep='\t').reset_index().iloc[:,1:]
        df_beta = df_beta[~df_beta.index.isin(remove)]
        if mode == 'pca':
            pca_loadings = pd.read_csv(
                os.path.join(spacia_res_path, 'model_input','sender_pc.csv'),index_col=0)
            df_beta = pd.DataFrame(
                np.matmul(
                    df_beta.values,
                    pca_loadings.values),
                columns = pca_loadings.columns
                )
        else:
            df_beta.columns = pathway_beta.Sender_pathway.unique()
        for gene in df_beta.columns:
            arr = []
            for i in range(n_chains):
                chain_df = df_beta[chain_size*i:chain_size*(i+1)][gene]
                arr += chain_df.sample(50, replace=False).tolist()
            arr = np.array(arr)
            mean_beta = df_beta[gene].mean()
            pval = stats.ttest_1samp(arr, 0)[1]
            if mode == 'pca':
                for key, val in zip(
                    pathway_beta.keys(),
                    [rg,gene,mean_beta,pval]
                    ):
                    pathway_beta[key].append(val)
            else:
                pathway_beta.loc[
                    (pathway_beta.index==rg) & (pathway_beta.Sender_pathway==gene),
                    'pval'
                    ] = pval
    if mode == 'pca':
        pathway_beta = pd.DataFrame(pathway_beta).set_index('RG')
    pathway_beta['pval_adj'] = p_adjust_bh(pathway_beta['pval'])
    return pathway_beta

def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]
#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Main function for running spacia, which evaluates interactions within \
            the context of cell neighborhoods, where the 'receiver' cells are the cells of interest, \
            and the cells from the neighborhood are referred to as 'sender' cells. The interactant \
            expressed in the receiver cells, through which the interactions are to be studied, are referred \
            to as 'Response', while the interactant expressed in the sender cells that potentially \
            influences the responder genes are called 'Signal'. The Goal of Spacia is to determine the \
            possibility and maginitude of how 'signal' gene(s) in 'sender' cells affect \
            'response' gene(s) in 'receiver' cells on a single-cell level."
)

    parser.add_argument(
        "counts",
        help="Path for gene expression data, spots by genes, must be normlized. \
            TXT format",
    )

    parser.add_argument(
        "spot_meta",
        help="Path for spot positional information, spots by feature(s). TXT format. \
            must have 'X', 'Y' columns for coordinates. 'cell_type' columns is needed \
            if running with '-rc' and '-sc' parameters. 'cell_type' refers to the group designation of cells, \
            e.g., type of cells. The user can specify which group of cells to use for spacia"
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
        "--receiver_features",
        "-rf",
        type=str,
        default=None,
        help="Input for the 'response' feature in the receiver cells. Can be: \
            1) a single gene as the 'response' feature. Works together with 'Corr_agg'; \
            2) multiple genes, separated by ',', each correspond to a different 'response' feature. \
            Works together with 'Corr_agg'; \
            3) multiple genes, sep by '|' used together as one 'response' feature. \
            'Corr_agg' if turned off in this mode; \
            4) a csv file, each row correponds to a 'response' feature. \
            The first column contains feature names, and the following columns are genes \
            that together form the feature. 'Corr_agg' if turned off in this mode;\
            5) 'pca', unsupervised mode, where the top ('num_comps') PCA embeddings are used as \
            'response' feature(s). 'Corr_agg' if turned off in this mode; \
            6) None, unsupervised mode, where the averaged expressions of gene clusters \
            in the 'receiver' cells are used as 'response' feature(s)."
    )

    parser.add_argument(
        "--sender_features",
        "-sf",
        type=str,
        default=None,
        help="Input for the 'signal' feature in the sender cells. Usage is similar to\
            '--receiver_features'",
    )

    parser.add_argument(
        "--corr_agg",
        "-ca",
        action="store_false",
        default=True,
        help="This is a toggle to turn off correlation based aggregation. Correlation aggregation \
            constructs 'response/signal' feature(s) based on expression of genes highly correlated \
            with the input gene. If it is turned off, the 'response/signal' feature(s) equal to the \
            expression of the input gene(s) ",
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
        "--corr_agg_method",
        "-cm",
        default='simple', # simple or weighted
        help="How the receiver gene expression will be aggregated, if 'simple', expression of top positively \
            correlated genes will be averaged. If 'weighted', expression of top genes ranked by absolute \
            correlation will be averaged using the correlation values as weights.",
    )
    
    parser.add_argument (
        "--num_comps",
        "-n_pc",
        help = "Number of components to be used in the pca mode.",
        default=20,
        type=int
    )
    
    parser.add_argument(
        "--response_exp_cutoff",
        "-rec",
        type=str,
        default=0.5,
        help="Quantile cutoff used to threshold the 'response' feature(s). If 'auto', a bimodel \
            distribution will be fitted to find the cutoff, and the cutoff will be calculated as \
            the average of the two peak means.",
    )
    
    parser.add_argument (
        "--pca_gene",
        "-pg",
        help = "Additional gene to be included in the pca mode. Note this will cause spacia \
            to use only the 3 top pcs. Experimental."
    )
    
    parser.add_argument(
        "--dist_cutoff",
        "-d",
        type=float,
        default=None,
        help="Distance cutoff for finding potential interacting cells.",
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

    parser.add_argument(
        "--output_path", "-o", type=str, default="spacia", help="Output path"
    )
    parser.add_argument(
        "--workers", type=int, default=min(4, os.cpu_count() or 1),
        help="Maximum concurrent R jobs; limited to the number of receiver features.",
    )
    parser.add_argument(
        "--cpp-cache-dir", default=None,
        help="Persistent compiled-model cache (default: .spacia_cache beside the code).",
    )
    parser.add_argument(
        "--zero-policy", choices=['keep', 'exclude_double_zero'], default='keep',
        help="Keep all spatial pairs (default), or remove pairs with zero raw transcripts in both named genes.",
    )
    parser.add_argument(
        "--detection-counts", default=None,
        help="Raw integer cell-by-gene count TSV for exclusion. Requires single-gene -sf/-rf; never uses scaled expression or binary responses to define zero.",
    )
    ######## Setting up ########
    # Debug param
    # args = parser.parse_args(
    #     [
    #         '/project/shared/xiao_wang/projects/cell2cell_inter/data/cosmx/results/nl_cpm.txt',
    #         '/project/shared/xiao_wang/projects/cell2cell_inter/data/cosmx/results/nl_cpm_metadata.txt',
    #         '-rc', 'Hep',
    #         '-sc', 'Inflammatory.macrophages',
    #         '-rf', 'nl_Hep_pathways.csv',
    #         '-sf', 'nl_Inflammatory.macrophages_pathways.csv',
    #         '-n', '25',
    #         '-b', '2',
    #         # '-nc', '20',
    #         # '-o', '/endosome/work/InternalMedicine/s190548/Spacia/test'
    #         ])
#%%
    ######## Setting up ########
    args = parser.parse_args()
    detection_counts, zero_config = prepare_zero_policy(args, parser)
    if args.cpp_cache_dir:
        os.environ["SPACIA_CACHE_DIR"] = os.path.abspath(args.cpp_cache_dir)
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
    response_exp_cutoff = args.response_exp_cutoff
    response_exp_cutoff = response_exp_cutoff if response_exp_cutoff == 'auto' else float(response_exp_cutoff)
    mcmc_params = args.mcmc_params
    corr_agg = args.corr_agg
    ntotal, nwarm, nthin, nchain = [int(x) for x in mcmc_params.split(",")]
    keep = args.keep_intermediate
    plot_mcmc = 'T' if args.plot_mcmc else 'F'
    ext = args.ext
    plot_debug = args.debug_plots
    corr_agg_method = args.corr_agg_method
    bag_size = args.bag_size
    nb = args.number_bags
    pca_gene = args.pca_gene
    n_pc = args.num_comps
    np.random.seed(0)

    # Checking inputs
    assert corr_agg_method in ['simple','weighted'], "'corr_agg_method' must be either 'simple' or 'weighted'!"  
    
    intermediate_folder = os.path.join(output_path, "model_input")
    if not os.path.exists(intermediate_folder):
        os.makedirs(intermediate_folder)
    dist_sender_fn = os.path.join(intermediate_folder, "dist_sender.json")
    metadata_fn = os.path.join(intermediate_folder, "metadata.txt")
    exp_sender_fn = os.path.join(intermediate_folder, "exp_sender.json")
    
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
    if args.zero_policy == 'exclude_double_zero' and any(
            gene not in cpm.columns for gene in [sender_features, receiver_features]):
        raise ValueError("Both selected detection genes must also exist in the expression matrix")

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
        
    # find candidate receiver and sender cells
    if dist_cutoff is None:
        dist_cutoff = calculate_neighbor_radius(
            spot_meta.iloc[:, :2], r_cells, s_cells, target_n_neighbors=n_neighbors, 
        )
        print(
            "Maximal distance for {} expected neighbors is {:.2f}".format(
                n_neighbors, dist_cutoff
            )
    )

    r2s_matrix = find_sender_candidates(
        r_cells, s_cells, spot_meta[["X", "Y"]], dist_cutoff
    )
    receiver_cell_for_cutoff = r2s_matrix.index.tolist()
    print('Limiting bags to those with at least {} sender cells'.format(bag_size))
    r2s_matrix = r2s_matrix[r2s_matrix.apply(len) >= bag_size]
    print('Number of bags: {}'.format(r2s_matrix.shape[0]))
    if r2s_matrix.shape[0] < 500 :
        # raise ValueError('Number of total bags is too small, job killed.')
        Warning('Number of total bags is too small.')
        pass
    elif r2s_matrix.shape[0]> nb:
        print('Subsample bags for Spacia.')
        r2s_matrix = r2s_matrix.sample(nb, replace=False)
    if r2s_matrix.empty:
        raise ValueError("No receiver neighborhoods satisfy the minimum sender count")
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
        n_pc,
        pca_gene
    )
    # If no receiver pathways are found, abort.
    if len(receiver_pathways.keys()) == 0:
        print('None of the genes in the provided receiver pathways are found in \
            the expression matrix, please modify the input and try again.')
        raise ValueError()
        
    print('Writing spacia_job.R inputs to the model_input folder.')
    # Calculate each receiver sender pair distances
    dist_r2s = r2s_matrix.to_frame().apply(
        lambda x: (cdist(
            # fixed issue in pandas that makes it object
            spot_meta.loc[[x.name], :"Y"].values.reshape(-1, 2),
            spot_meta.loc[x.iloc[0], :"Y"].values.reshape(-1, 2),
        )[0]/dist_cutoff).round(5), # normalize distance to 0-1
        axis=1,
    )
    
    sender_dist_dict = {}
    for i in dist_r2s.index:
        sender_dist_dict[i] = dist_r2s[i].tolist()

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
        
    # Fit features on the original cohort, then apply the raw-count pair mask.
    # The receiver cutoff population above is also left unchanged.
    distance_normalizer = None
    if args.zero_policy == 'exclude_double_zero':
        original_r2s = r2s_matrix
        distance_normalizer = max(max(d) for d in sender_dist_dict.values())
        r2s_matrix, pair_audit, receiver_audit = filter_sender_pairs(
            original_r2s, args.zero_policy, detection_counts,
            args.sender_features, args.receiver_features, bag_size)
        pair_audit.to_csv(os.path.join(output_path, 'zero_policy_pairs.csv'), index=False)
        receiver_audit.to_csv(os.path.join(output_path, 'zero_policy_receivers.csv'), index=False)
        zero_summary = dict(zero_config,
            receivers_before=len(original_r2s), receivers_retained=len(r2s_matrix),
            receivers_removed=len(original_r2s)-len(r2s_matrix),
            pairs_before=len(pair_audit),
            double_zero_pairs_removed=int(pair_audit.double_zero.sum()),
            pairs_removed_small_neighborhood=int(pair_audit.action.eq('excluded_small_neighborhood').sum()),
            pairs_retained=int(pair_audit.action.eq('retained').sum()),
            distance_normalizer=distance_normalizer)
        with open(os.path.join(output_path, 'zero_policy_summary.json'), 'w') as handle:
            json.dump(zero_summary, handle, indent=2)
        if r2s_matrix.empty:
            raise ValueError("Double-zero exclusion leaves no eligible receiver neighborhoods; see zero_policy_* audits")
        if not np.isfinite(distance_normalizer) or distance_normalizer <= 0:
            raise ValueError("Cannot normalize distances: the original neighborhoods have no positive distance")
        filtered_distances = {}
        for receiver, senders in r2s_matrix.items():
            keep_senders = set(senders)
            filtered_distances[receiver] = [d for sender, d in zip(
                original_r2s[receiver], sender_dist_dict[receiver]) if sender in keep_senders]
        sender_dist_dict = filtered_distances
        receiver_candidates = r2s_matrix.index.tolist()
        retained_senders = {s for senders in r2s_matrix for s in senders}
        sender_candidates = [s for s in sender_candidates if s in retained_senders]

    # contruct and save metadata
    meta_data = spot_meta.loc[receiver_candidates, :"Y"]
    meta_data["Sender_cells"] = r2s_matrix.loc[receiver_candidates].apply(",".join)
    meta_data_senders = spot_meta.loc[sender_candidates, :"Y"]
    meta_data = pd.concat([meta_data, meta_data_senders])

    sender_exp = (
        r2s_matrix.to_frame()
        .apply(lambda x: sender_pathway_exp.loc[x[0],].values.round(3).tolist(), axis=1)
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
        if job_finished and args.zero_policy == 'keep':
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
        # print(response_exp_cutoff)
        rf_to_drop = []
        if response_exp_cutoff == 'auto':
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
            cutoff = receiver_exp.quantile(response_exp_cutoff)
    
        if plot_debug:
            receiver_exp.hist(bins=20,density=True)
            plt.plot((cutoff,cutoff), (0,2))
            plt.savefig(
                os.path.join(intermediate_folder, job_id + "_exp_receiver_dist.pdf"))
            plt.close()
            
        receiver_exp = receiver_exp > cutoff
        receiver_exp = receiver_exp + 0
        if args.zero_policy == 'exclude_double_zero':
            receiver_audit['response_binary'] = receiver_audit.receiver.map(receiver_exp)
            receiver_audit.to_csv(os.path.join(output_path, 'zero_policy_receivers.csv'), index=False)
            before = receiver_exp.loc[original_r2s.index]
            after = receiver_exp.loc[receiver_candidates]
            zero_summary.update(response_cutoff=float(cutoff),
                response_0_before=int(before.eq(0).sum()), response_1_before=int(before.eq(1).sum()),
                response_0_retained=int(after.eq(0).sum()), response_1_retained=int(after.eq(1).sum()))
            with open(os.path.join(output_path, 'zero_policy_summary.json'), 'w') as handle:
                json.dump(zero_summary, handle, indent=2)
        receiver_exp = receiver_exp[receiver_candidates]
        receiver_exp.to_csv(exp_receiver_fn, header=None, index=None)

        if job_finished:
            print(job_id + ' is already finished; refreshed zero-policy audits without refitting.')
            continue
        spacia_output_path = os.path.join(output_path, job_id)
        if not os.path.exists(spacia_output_path):
            os.makedirs(spacia_output_path)
        spacia_jobs.append(
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
                ] + (["1", str(distance_normalizer)] if distance_normalizer is not None else [])
        )
    
    with open(os.path.join(output_path, 'spacia_r.log'), 'w') as f:
        f.write('\n'.join(shlex.join(job) for job in spacia_jobs)) # Save the actual jobs for debug purpose
        
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

        with open(os.path.join(intermediate_folder, fn), "w") as f:
            f.write(format_json(pathway_dict))
            
    # Writing spacia R job inputs common for each receiver pathways
    # job metadata
    meta_data.to_csv(metadata_fn, sep='\t')
    
    # sender distance and expression json (list of lists)
    with open(dist_sender_fn, "w") as f:
        f.write(format_json(sender_dist_dict))
        
    # with open(exp_sender_fn, "w") as f:
    #     f.write(format_json(sender_exp))
    with open(exp_sender_fn, "w") as f:
        f.write(format_json(sender_exp))
    
    ######## Proceed with spacia_job.R ########
    # Run all spacia R jobs
    print('Running spacia_R MCMC MIL models.')
    if args.workers < 1:
        raise ValueError("--workers must be positive")
    if spacia_jobs:
        # Each bounded worker processes a batch in the same R session. Threads
        # only supervise subprocesses; they do not duplicate Python/NumPy state.
        n_workers = min(args.workers, len(spacia_jobs))
        batch_commands = []
        for worker in range(n_workers):
            manifest = os.path.join(intermediate_folder, f"worker_{worker}.json")
            with open(manifest, "w") as f:
                json.dump([job[2:] for job in spacia_jobs[worker::n_workers]], f)
            batch_commands.append(["Rscript", os.path.join(spacia_path, "spacia_batch.R"), manifest])
        if n_workers == 1:
            spacia_worker(batch_commands[0])
        else:
            with ThreadPoolExecutor(max_workers=n_workers) as pool:
                list(pool.map(spacia_worker, batch_commands))

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
        res_beta.Sender_pathway = sender_pathways_names
        pathways = pd.concat([pathways, res_beta])

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
        interactions = pd.concat([interactions, _interactions])

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
        b_plus_fdr = pd.concat([b_plus_fdr, fdr])
        
    # update pathway_betas
    c_l = int((ntotal-nwarm)/nthin)
    agg_mode = 'pca' if sender_features == 'pca' else 'gene'
    pathways = process_beta(pathways.copy(), output_path, c_l, nchain,agg_mode)
    pathways.to_csv(os.path.join(output_path, "Pathway_betas.csv"))
    
    interactions.to_csv(os.path.join(output_path, "Interactions.csv"))
    # calculate p values for b
    b_plus_fdr = process_b(b_plus_fdr.copy(), output_path, c_l, nchain)
    b_plus_fdr.to_csv(os.path.join(output_path, "B_and_FDR.csv"))
    
    # Remove model_input files
    if not keep:
        os.system("rm -rf {}".format(intermediate_folder))
