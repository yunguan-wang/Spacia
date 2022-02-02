"""
# TODO: Fillme

Dependency
----------
R >= 4.0.2
    Rcpp

Python dependencies will be handeled automatically by installing the scia package.

"""

from dataclasses import replace
from itertools import count
import os
import sys
import argparse
import logging
from multiprocessing import Pool
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
import subprocess


def scia_worker(cmd):
    """
    # TODO: Fillme
    """
    # NOTE: fill code
    pass
    # os.system('rm -r output_fn') # Get rid of output.
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
        description="Main function for running SCIA-MIL",
    )

    parser.add_argument(
        "input_path", type=str, 
        help="Input folder containing all necessary files."
    )

    parser.add_argument(
        "output_path", type=str, help="Output path")

    parser.add_argument(
        "--counts",
        "-su",
        help="Gene expression data, spots by genes. TXT format",
    )

    parser.add_argument(
        "--spot_meta",
        "-s",
        help="spot positional information, spots by features. TXT format",
    )

    parser.add_argument(
        "--n_neighbors",
        "-n",
        type = float,
        default = 10,
        help= "Expected number of nearest neighbors for most of the cells. The \
        exact number could vary from cell to cell.",
    )
    # parser.add_argument(
    #     "--scia_R",
    #     "-r",
    #     default=0.08,
    #     type=float,
    #     help="Spot neighborhood radius ratio, 0-1, \
    #         radius=R*min(xmax-xmin,ymax-ymin).",
    # )

# Setting up
####
    input_path = '/project/shared/xiao_wang/projects/cell2cell_inter/data'
    output_path = '/project/shared/xiao_wang/projects/cell2cell_inter/data'
    pathway_lib = 'pathways.txt'
    n_neighbors = 10
    counts = 'Counts.txt'
    spot_meta = 'Spot_metadata.txt'
####
    # args = parser.parse_args()
    # input_path = args.input_path
    # output_path = args.output_path
    # counts = args.counts
    # spot_meta = args.spot_meta

    counts = pd.read_csv(
        os.path.join(input_path, counts), index_col=0, sep='\t'
    )
    # TODO: add count normalizations.
    spot_meta = pd.read_csv(
        os.path.join(input_path, spot_meta), index_col=0, sep='\t'
    )
    pathway_lib = pd.read_csv(
        os.path.join(input_path, pathway_lib), index_col=0, sep='\t'
    )
    # getting script path for supporting codes.
    scia_path = os.path.abspath(__file__)
    scia_path = "/".join(scia_path.split("/")[:-1]) + '/scia'
    # scia_script = os.path.join(scia_path, "denoise.R")

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Setting up logs
    log_fn = os.path.join(output_path, "scia_log.txt")
    if os.path.exists(log_fn):
        os.remove(log_fn)
    logging.basicConfig(
        filename=log_fn,
        format="%(asctime)s,%(levelname)s:::%(message)s",
        datefmt="%H:%M:%S",
        level="INFO",
    )
    # redirects stdout and stderr to logger
    stdout_logger = logging.getLogger("STDOUT")
    sl = StreamToLogger(stdout_logger, logging.INFO)
    sys.stdout = sl
    stderr_logger = logging.getLogger("STDERR")
    sl = StreamToLogger(stderr_logger, logging.ERROR)
    sys.stderr = sl

####
    neighbor_r = calculate_neighbor_radius(spot_meta)
    receivers = np.random.choice(counts.index, 100, replace=False)
    dist_matrix = cdist(spot_meta,spot_meta)
    np.fill_diagonal(dist_matrix, np.inf)
    dist_matrix = pd.DataFrame(
        dist_matrix,
        index = spot_meta.index, columns = spot_meta.index)
    senders = dist_matrix.loc[receivers]<=neighbor_r
    senders = senders[senders.any(axis=1)]
    senders = senders.apply(lambda x: x[x].index.tolist(), axis=1)
    # Define the matrix shared by all child R processes
    sender_cells = np.unique(senders.sum())

    receiver_genes = pathway_lib.index.unique()
    for receiver_gene in receiver_genes:
        sender_genes = pathway_lib.loc[receiver_gene].iloc[:,0].unique()
        sender_cts = counts.loc[sender_cells, sender_genes].round(2)
        sender_dist = dist_matrix.loc[receivers, sender_cells]
        recivers_mtx = pd.DataFrame(
            senders.apply(lambda x: ','.join(x))).rename_axis(None)
        recivers_mtx.columns = ['senders']
        exp_receiver = counts.loc[receivers, receiver_gene]
        # NOTE: place holder rule for dichtomizing receiver gene exp
        exp_receiver = 0 + (exp_receiver > exp_receiver.median())
        recivers_mtx['exp_receiver'] = exp_receiver
        # Write out R wrapper data
        for df, f_name in zip(
            [sender_cts, sender_dist, recivers_mtx],
            ['job_cts.txt', 'job_dist.txt', 'job_receiver.txt']):
            df.to_csv(
                os.path.join(output_path, f_name), sep='\t')
        # TODO: R job submitting, results collection and clean up

