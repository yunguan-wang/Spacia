"""
# TODO: Fillme

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
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import robust_scale


def spacia_worker(cmd):
    """
    worker function for multiprocessing.
    """
    os.system(cmd)
    # remove temp job input files
    # os.system('rm -f {}'.format(' '.join(spacia_job_inputs)))
    return


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


def find_sender_candidates(receivers, senders, locations, dist_cutoff=30):
    """
    receivers, senders: list of spot ids.
    locations: pd.DataFrame of X, Y locations
    """
    pip = pd.Series(dtype=object)
    for r in receivers:
        dist_to_senders = cdist(
            locations.loc[r].values.reshape(-1, 2),
            locations.loc[senders].values.reshape(-1, 2),
        )[0]
        crit = dist_to_senders <= dist_cutoff
        if crit.sum() == 0:
            continue
        pip[r] = senders[crit].tolist()
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
                keep_genes = ~crit
    filtered_counts = counts.iloc[keep_cells.values, keep_genes.values]
    filtered_cpm = filtered_counts.apply(lambda x: 1e4 * x / x.sum(), axis=1)
    return filtered_cpm

def get_corr_agg_genes(corr_agg, cpm, cells, g, top_corr_genes):
    if corr_agg:
        print('Constructing pathway using correlation aggregation')
        corr = cdist(
            cpm.loc[cells, g].values.reshape(1, -1),
            cpm.loc[cells].T,
            metric="correlation",
        )[0]
        pathway_genes = cpm.columns[np.argsort(corr)[:top_corr_genes]]
        pathway_name = g + "_correlated genes"
    else:
        logging.warning("Correlation aggregation is turned off and \
            this pathway has only one gene. This is not recommended.")
        pathway_name = g
    return pathway_genes, pathway_name

def contruct_pathways(
    cpm, receiver_candidates, sender_candidates, receiver_features, sender_features,
):
    receiver_pathways = {}
    sender_pathways = {}
    for pathway_dict, pathway_features, pathway_type in zip(
        [receiver_pathways, sender_pathways],
        [receiver_features, sender_features],
        ["Receiver", "Sender"],
    ):
        cells  = receiver_candidates if pathway_type == "Receiver" else sender_candidates
        if pathway_features is None:
            print(
                "{} features is not provides, use gene modules as sender pathways.".format(pathway_type)
            )
            # calculate clusters
            pathway_exp = cpm.loc[
                cells,
            ]
            pathway_exp = pathway_exp.T[pathway_exp.std() > 0].T
            # add an expression level cutoff to reduce the number of genes
            top_expressed_genes = pathway_exp.mean().sort_values()[
                -int(pathway_exp.shape[1] / 4) :
            ]
            pathway_exp = pathway_exp[top_expressed_genes.index]
            gene_clusters = AgglomerativeClustering(
                50, affinity="correlation", linkage="complete"
            ).fit_predict(robust_scale(pathway_exp.T))
            # construct sender_pathway
            print("Cosntruct 50 {} pathways from gene modules".format(pathway_type))
            for cluster in np.unique(gene_clusters):
                gene_mask = gene_clusters == cluster
                pathway_dict["module_" + str(cluster + 1)] = pathway_exp.columns[
                    gene_mask
                ].tolist()
            continue
        if pathway_features[-4:] == ".csv":
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
                            corr_agg, cpm, cells, g, top_corr_genes)
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
                    corr_agg, cpm, cells, g, top_corr_genes)
                pathway_dict[pathway_name] = pathway_genes.tolist()
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


#%%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Main function for running spacia. This are running modes and spacia will \
            automatically decide which one to use based on the 'receiver/sender features' and \
            the 'corr_agg' toggle. Below are the four modes: \n\t\
            (1) No aggregation: The user either provides one or several single genes, \
            for sending/receiving pathways. These will be used in the analysis directly; \n\t\
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
        help="Path for gene expression data, spots by genes, normalized. \
            TXT format",
    )

    parser.add_argument(
        "spot_meta",
        help="Path for spot positional information, spots by features. TXT format. \
            must have 'X', 'Y' columns for coordinates. 'cluster' columns is needed \
            if running with '-rc' and '-sc' parameters. 'cluster' refers to the group designation of cells, \
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
            5) if nothing is passed, spacia will run in mode (4).",
    )

    parser.add_argument(
        "--sender_features",
        "-sf",
        type=str,
        default=None,
        help="Sender gene(s). Can be: \
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
            5) if nothing is passed, spacia will run in mode (4).",
    )

    parser.add_argument(
        "--output_path", "-o", type=str, default="spacia", help="Output path"
    )

    parser.add_argument(
        "--dist_cutoff",
        "-d",
        type=str,
        default=None,
        help="Distance cutoff for finding potential interacting cells.",
    )

    parser.add_argument(
        "--receiver_exp_cutoff",
        "-rec",
        type=float,
        default=0.5,
        help="Quantile cutoff used to threshold receiver expression.",
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
        default="50000,20000,100,1",
        help="MCMC parameters, four values packed here are {ntotal,nwarm,nthin,nchain}",
    )


    parser.add_argument(
        "--receiver_cluster",
        "-rc",
        type=str,
        default=None,
        help="Name of receiver cluster, must be in spot_metadata.",
    )

    parser.add_argument(
        "--sender_cluster",
        "-sc",
        type=str,
        default=None,
        help="Name of sender cluster, must be in spot_metadata.",
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
        help="Whether the intermediate folder should be deleted.",
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

    ######## Setting up ########
    # Debug params
    # args = parser.parse_args(
    #     [
    #     '/project/shared/xiao_wang/projects/cell2cell_inter/data/Counts.txt',
    #     '/project/shared/xiao_wang/projects/cell2cell_inter/data/Spot_metadata.txt',
    #     'FGFR1',
    #     '-o', '/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/spacia_py_test',
    #     '-rc', 'C_1',
    #     '-sc', 'C_2',
    #     '-cf', '/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/spacia_py_test/input/input_cells.csv',
    #     ]
    # )

    ######## Setting up ########
    args = parser.parse_args()
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
    mcmc_params = args.mcmc_params
    corr_agg = args.corr_agg
    ntotal, nwarm, nthin, nchain = mcmc_params.split(",")
    keep = args.keep_intermediate

    intermediate_folder = os.path.join(output_path, "intermediate")
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

    # getting script path for supporting codes.
    spacia_path = os.path.abspath(__file__)
    wrapup_script = spacia_path.replace("spacia.py", "wrapup_spacia_results.py")
    spacia_path = "/".join(spacia_path.split("/")[:-2]) + "/spacia"
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
    cpm = preprocessing_counts(counts)
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
    if (sender_cluster is not None) & (sender_cluster is not None):
        r_cells = spot_meta[spot_meta.cluster == receiver_cluster].index
        s_cells = spot_meta[spot_meta.cluster == sender_cluster].index
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
    sender_candidates = list(set(r2s_matrix.sum()))
    receiver_candidates = r2s_matrix.index.tolist()

    ######## Preparing spacia_job.R inputs ########
    # Contruct sender and receiver pathways
    receiver_pathways, sender_pathways = contruct_pathways(
        cpm, receiver_candidates, sender_candidates, receiver_features, sender_features
    )
    # If no receiver pathways are found, abort.
    if len(receiver_pathways.keys()) == 0:
        print('None of the genes in the provided receiver pathways are found in \
            the expression matrix, please modify the input and try again.')
        raise ValueError()
        
    for pathway_dict, fn in zip(
        [receiver_pathways, sender_pathways],
        ["receiver_pathways.json", "sender_pathways.json"],
    ):
        with open(os.path.join(intermediate_folder, fn), "w") as fp:
            json.dump(pathway_dict, fp) # Save receiver and sender pathways

    print('Writing spacia_job.R inputs to the intermediate folder.')
    dist_sender_fn = os.path.join(intermediate_folder, "dist_sender.json")
    metadata_fn = os.path.join(intermediate_folder, "metadata.csv")
    exp_sender_fn = os.path.join(intermediate_folder, "exp_sender.json")
    # Calculate each receiver sender pair distances
    dist_r2s = r2s_matrix.to_frame().apply(
        lambda x: cdist(
            spot_meta.loc[x.name, :"Y"].values.reshape(-1, 2),
            spot_meta.loc[x[0], :"Y"].values.reshape(-1, 2),
        )[0].round(2),
        axis=1,
    )
    sender_dist_dict = {}
    for i in dist_r2s.index:
        sender_dist_dict[i] = dist_r2s[i].tolist()
    with open(dist_sender_fn, "w") as fp:
        json.dump(sender_dist_dict, fp)

    # contruct and save metadata
    meta_data = spot_meta.loc[receiver_candidates, :"Y"]
    meta_data["Sender_cells"] = r2s_matrix.loc[receiver_candidates].apply(",".join)
    meta_data.to_csv(metadata_fn, sep="\t")

    # contruct and save sender exp
    sender_pathway_exp = pd.DataFrame(
        index=sender_candidates, columns=sender_pathways.keys()
    )
    for key in sender_pathway_exp.columns:
        sender_pathway_exp[key] = (
            cpm.loc[sender_candidates, sender_pathways[key]].mean(axis=1).values
        )
    sender_exp = (
        r2s_matrix.to_frame()
        .apply(lambda x: sender_pathway_exp.loc[x[0],].values.round(2).tolist(), axis=1)
        .to_dict()
    )
    with open(exp_sender_fn, "w") as fp:
        json.dump(sender_exp, fp)

    ######## Write spacia_job.R jobs ########
    # construct receiver expression and the job commands
    spacia_jobs = []
    spacia_job_folders = []
    for rp in receiver_pathways.keys():
        job_id = rp
        spacia_job_folders.append(os.path.join(output_path, job_id))
        exp_receiver_fn = os.path.join(
            intermediate_folder, job_id + "_exp_receiver.csv"
        )
        # Getting receiver exp
        rp_genes = receiver_pathways[rp]
        receiver_exp = cpm.loc[receiver_candidates, rp_genes].mean(axis=1)
        receiver_exp = receiver_exp > receiver_exp.quantile(receiver_exp_cutoff)
        receiver_exp = receiver_exp + 0
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
                ]
            )
        )
    
    with open(os.path.join(output_path, 'spacia_r.log'), 'w') as f:
        f.write('\n'.join(spacia_jobs)) # Save the actual jobs for debug purpose

    ######## Proceed with spacia_job.R ########
    # Run all spacia R jobs
    print('Running spacia_R MCMC MIL models.')
    with Pool(16) as p:
        _ = p.map(spacia_worker, spacia_jobs)
    
    # Collect all results
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
        res_beta = pd.read_csv(
            os.path.join(fd, job_id + "_beta.txt"), sep="\t"
        ).mean()
        res_beta = res_beta.reset_index()
        res_beta.index = [job_id] * res_beta.shape[0]
        res_beta.columns = ["Sender_pathway", "Beta"]
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
    
    # Remove intermediate files
    if not keep:
        os.system("rm -rf {}".format(intermediate_folder))