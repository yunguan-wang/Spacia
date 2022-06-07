"""
# TODO: Fillme

Dependency
----------
R >= 4.0.2
    Rcpp
if running in biohpc, need to do module purge; module load shared first.
Python dependencies will be handeled automatically by installing the spacia package.

Test:
cd /project/shared/xiao_wang/projects/cell2cell_inter/data
python ../code/spacia/spacia.py mini_Counts.txt Spot_metadata.txt C_1 C_2 pathways.txt
"""

from itertools import count
import os
import sys
import argparse
import logging
from multiprocessing import Pool
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

def spacia_worker(spacia_job_params):
    """
    worker function for multiprocessing.
    """
    spacia_job_inputs = spacia_job_params.split(' ')[2:5]
    os.system('Rscript ' + spacia_job_params)
    # remove temp job input files
    # os.system('rm -f {}'.format(' '.join(spacia_job_inputs)))
    return

def calculate_neighbor_radius(
    spot_meta, sample_size=1000, target_n_neighbors = 10, margin = 10,
    r_min = 1, r_max = 1000):
    r_samples = spot_meta.loc[
        np.random.choice(spot_meta.index,sample_size, False)].copy().loc[:,['X','Y']]
    r_sample_dist = cdist(r_samples, spot_meta.loc[:,['X','Y']])
    n_steps = 1
    while n_steps <= 100:
        r_next = (r_max + r_min)/2
        n_neighbors = np.median(np.sum((r_sample_dist <= r_next),axis=1))
        if n_neighbors == target_n_neighbors:
            break
        elif (r_max - r_min)/2 < margin:
            break 
        n_steps += 1
        nn_r_min = np.median(np.sum((r_sample_dist <= r_min),axis=1))
        nn_r_max = np.median(np.sum((r_sample_dist <= r_max),axis=1))
        if (n_neighbors>target_n_neighbors) == (nn_r_max > target_n_neighbors):
            r_max = r_next
        elif (n_neighbors>target_n_neighbors) == (nn_r_min> target_n_neighbors):
            r_min = r_next
    return r_next

def find_sender_candidates(receivers, senders, locations, dist_cutoff = 30):
    '''
    receivers, senders: list of spot ids.
    locations: pd.DataFrame of X, Y locations
    '''
    pip = pd.Series(dtype = object)
    for r in receivers:
        dist_to_senders = cdist(
            locations.loc[r].values.reshape(-1,2), 
            locations.loc[senders].values.reshape(-1,2))[0]
        crit = dist_to_senders <= dist_cutoff
        if crit.sum() == 0:
            continue
        pip[r] = senders[crit].tolist()
    return pip

def preprocessing_counts(
    counts, ntotal_cutoff = 100, n_genes_cutoff = 20, n_cells_cutoff = 10):
    '''
    Simple QC based on total counts, num_genes expressed in cell and num of cells a gene is expressed.
    '''
    counts = counts.T.groupby(counts.columns).mean().T
    # add the total counts per cell as observations-annotation to counts
    n_total = counts.sum(axis=1)
    n_genes = np.sum(counts>0,axis=1)
    n_cells = np.sum(counts>0,axis=0)
    keep_cells = np.zeros(shape=[counts.shape[0]], dtype=bool)
    keep_cells[:] = True
    for i, df, cutoff, qc_type in zip(
        [0, 0, 1],
        [n_total, n_genes, n_cells],
        [ntotal_cutoff, n_genes_cutoff, n_cells_cutoff],
        ['Total counts', 'Total genes', 'Number of cells expressed in']):
        crit = df < cutoff
        num_bad = crit.sum()
        if num_bad > 0:
            print('{} {} are dropped because {} is less than {}.'.format(
                num_bad, 'Cells' if i==0 else 'Genes', qc_type, cutoff
            ))
            if i == 0:
                keep_cells = np.logical_and(keep_cells, ~crit)
            else:
                keep_genes = ~crit
    filtered_counts = counts.iloc[keep_cells.values, keep_genes.values]
    filtered_cpm = filtered_counts.apply(lambda x: 1e4*x/x.sum(), axis=1)
    return filtered_cpm

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



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Main function for running spacia",
    )

    parser.add_argument(
        "counts",
        help="Path for gene expression data, spots by genes, normalized. \
            TXT format",
    )

    parser.add_argument(
        "spot_meta",
        help="Path for spot positional information, spots by features. TXT format. \
            must have 'X', 'Y', and 'cluster' columns",
    )

    parser.add_argument(
        "receiver_cluster", type=str,
        help="Name of receiver cluster, must be in spot_metadata."
    )    # TODO: Add options to take a txt file of cell ids instead of clusters.

    parser.add_argument(
        "sender_cluster", type=str,
        help="Name of sender cluster, must be in spot_metadata."
    )    # TODO: Add options to take a txt file of cell ids instead of clusters.

    parser.add_argument(
        "receiver_features", 
        type = str,
        help = "Receiver gene(s). Can be 1) a single gene; 2) multiple genes, separated by ',' \
            each is treated as the seed of receiver pathway; 3) multiple genes, sep by '|' \
            used as one pathway including provided genes; or 4) a csv file each row for a \
            receiver pathway and columns are genes ."
    )   # TODO: To be implemented

    parser.add_argument(
        "--sender_features",
        "-sf",
        type = str,
        defaule = None,
        help = "Sender gene(s). Can be 1) a single gene; 2) multiple genes, separated by ',' \
            each is treated as the seed of receiver pathway; 3) multiple genes, sep by '|' \
            used as one pathway including provided genes; or 4) a csv file each row for a \
            sender pathway and columns are genes, must match the csv file of the receiver \
            features input option 4)." 
    )   # TODO: To be implemented
    
    parser.add_argument(
        "--output_path", 
        "-o",
        type=str,
        default = 'spacia',
        help="Output path"
    )

    parser.add_argument(
        "--dist_cutoff", 
        "-d", 
        type=str,
        default=None,
        help="Name of sender cluster, must be in spot_metadata."
    )

    parser.add_argument(
        "--num_corr_genes", 
        "-nc", 
        type = int,
        default = 100,
        help = "Number of correlated gene to use in calculating receiver pathway expression."
    )

    parser.add_argument(
        "--n_neighbors",
        "-n",
        type = float,
        default = 10,
        help= "Expected number of nearest neighbors for most of the cells. The \
        exact number could vary from cell to cell.",
    )

    parser.add_argument(
        "--mcmc_params",
        "-m",
        type = str,
        default = '8000,4000,10,1',
        help= "MCMC parameters, four values packed here are {ntotal,nwarm,nthin,nchain}"
    )


######## Setting up ########

    # Debug params
    args = parser.parse_args(
        [
        '/project/shared/xiao_wang/projects/cell2cell_inter/data/Counts.txt',
        '/project/shared/xiao_wang/projects/cell2cell_inter/data/Spot_metadata.txt',
        'C_1',
        'C_2',
        'FGFR1',
        '/project/shared/xiao_wang/projects/cell2cell_inter/data/pathways.txt']
        )

    args = parser.parse_args()
    counts = args.counts
    spot_meta = args.spot_meta
    receiver_cluster = args.receiver_cluster
    sender_cluster = args.sender_cluster
    pathway_lib = args.pathway_lib
    output_path = args.output_path
    n_neighbors = args.n_neighbors
    receiver_features = args.receiver_features
    top_corr_genes = args.num_corr_genes
    dist_cutoff = args.dist_cutoff
    mcmc_params = args.mcmc_params.replace(',', ' ')

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
    
    counts = pd.read_csv(counts, index_col=0, sep='\t')
    spot_meta = pd.read_csv(spot_meta, index_col=0, sep='\t')
    cpm = preprocessing_counts(counts)
    cpm, spot_meta = cpm.align(spot_meta, join = 'inner', axis=0)
    # pathway_lib = pd.read_csv(pathway_lib, index_col=0, sep='\t')
    # getting script path for supporting codes.
    spacia_path = os.path.abspath(__file__)
    spacia_path = "/".join(spacia_path.split("/")[:-1]) + '/spacia'
    spacia_script = os.path.join(spacia_path, "spacia_job.R")

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # redirects stdout and stderr to logger
    stdout_logger = logging.getLogger("STDOUT")
    sl = StreamToLogger(stdout_logger, logging.INFO)
    sys.stdout = sl
    stderr_logger = logging.getLogger("STDERR")
    sl = StreamToLogger(stderr_logger, logging.ERROR)
    sys.stderr = sl

######## Preparing spacia jobs ########
    if dist_cutoff is None:
        dist_cutoff = calculate_neighbor_radius(
            spot_meta.iloc[:,:2], target_n_neighbors=n_neighbors)
    print(
        'Maximal distance for {} expected neighbors is {:.2f}'.format(
            n_neighbors, dist_cutoff
        ))
    receiver_candidates = spot_meta[spot_meta.cluster == receiver_cluster].index
    senders_candidates = spot_meta[spot_meta.cluster == sender_cluster].index
    dist_matrix = cdist(spot_meta.iloc[:,:2],spot_meta.iloc[:,:2])
    r2s_matrix = find_sender_candidates(
        receiver_candidates,
        senders_candidates,
        spot_meta[['X','Y']],
        dist_cutoff)

    # Getting receiver expression
    if receiver_features[-4:] == '.csv':
        pass
    elif '|' in receiver_features:
        pass
    elif ',' in receiver_features:
        pass
    else:
        receiver_gene = receiver_features
        corr = cdist(
            cpm.loc[r2s_matrix.index, receiver_gene].values.reshape(1,-1),
            cpm.loc[r2s_matrix.index].T,
            metric = 'correlation'
        )[0]
        top_genes = cpm.columns[np.argsort(corr)[:100]]
    #!TODO : continue here
    # TODO:


    # Define the matrix shared by all child R processes
    sender_cells = np.unique(senders.sum())
    receiver_genes = pathway_lib.index.unique()
    print('Total number of receiver cells: {}'.format(senders.shape[0]))
    print('Total number of sender cells: {}'.format(len(sender_cells)))
    # generate all spacia jobs
    spacia_r_jobs_list = []
    for _, receiver_gene in enumerate(receiver_genes):
        print('Processing {}'.format(receiver_gene))
        sender_genes = pathway_lib.loc[receiver_gene].iloc[:,0].unique()
        sender_cts = counts.loc[sender_cells, sender_genes].round(2)
        sender_dist = dist_matrix.loc[receivers, sender_cells]
        sender_dist = (sender_dist / dist_cutoff).round(2)
        recivers_mtx = pd.DataFrame(
            senders.apply(lambda x: ','.join(x))).rename_axis(None)
        recivers_mtx.columns = ['senders']
        exp_receiver = counts.loc[receivers, receiver_gene]
        # NOTE: place holder rule for dichtomizing receiver gene exp
        exp_receiver = 0 + (exp_receiver > exp_receiver.median())
        recivers_mtx['exp_receiver'] = exp_receiver
        job_id = '_'.join(['spacia_job', receiver_gene])
        # Write out R wrapper data
        spacia_job_inputs = []
        for df, f_name in zip(
            [sender_cts, sender_dist, recivers_mtx],
            ['_cts.txt', '_dist.txt', '_receiver.txt']):
            df.to_csv(
                os.path.join(output_path, job_id+f_name), sep='\t')
            spacia_job_inputs.append(
                os.path.abspath(output_path) + '/' + job_id+f_name
            )
        # construct spacia.R jobs
        spacia_job_params = ' '.join([
            spacia_script,
            spacia_script.replace('spacia_job.R',''),
            ' '.join(spacia_job_inputs),
            job_id,
            mcmc_params,
            output_path + '/'
            ]
            )
        spacia_r_jobs_list.append(spacia_job_params)
    # run all jobs using 16 processes
    with Pool(16) as p:
        _ = p.map(spacia_worker, spacia_r_jobs_list)
