'''
Requires: R/4.0.2, gcc
'''
#%%
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from collections import defaultdict
import json
import os
import subprocess
from multiprocessing import Pool
import sys

def simulate_grid(grid_size, num_dots, spacing=15):
    locations = np.random.uniform(0, grid_size, size=(1,2))
    iter = 0
    while locations.shape[0]<num_dots:
        iter += 1
        if iter >= 100000:
            print(
                'Max iter reached, final dataset size is {}'.format(
                    locations.shape[0]
                    )
                )
            break
        tmp = np.random.uniform(0, grid_size, size=(1,2))
        dist = cdist(tmp, locations)[0]
        if (dist<=spacing).any():
            continue
        else:
            locations = np.append(locations, tmp, axis=0)
    return locations

def estimate_radius(v, center, r_min, r_max, shape = 'blob'):
    shape = {'blob':'blob', 'ring':'ring'}[shape]
    if shape == 'blob':
        r = np.sqrt(sum((v-center)**2))
        return r < r_max
    elif shape == 'ring':
        r = np.sqrt(sum((v-center)**2))
        return (r > r_min) & (r < r_max)

def simulate_pip(receivers, senders, locations, dist_cutoff = 30):
    pip = []
    senders = np.array(senders)
    for r in receivers:
        dist_to_senders = cdist(locations[[r]], locations[senders])[0]
        crit = dist_to_senders <= dist_cutoff
        pip.append(senders[crit])
    return np.array(pip,dtype=object)

def simulate_senders_receivers(
    num_dots, locations,
    centers = np.array([[292, 799],[894, 149],[855, 774],[611, 391],[297, 366]]),
    radius_min = 75, radius_max = 150, ring_r = 50):
    np.random.seed(0)
    radius = np.random.randint(radius_min, radius_max, size = 5)
    senders = []
    receivers = []
    for c, R in zip(centers,radius):
        _sender = [
            i for i in range(num_dots) if estimate_radius(
                locations[i], c, None, R)
                ]
        _receiver = [
            i for i in range(num_dots) if estimate_radius(
                locations[i], c, R,R+ring_r,shape = 'ring')
                ]
        senders += _sender
        receivers += _receiver
    senders = list(set(senders))
    receivers = list(set(receivers))
    return senders, receivers

def simulate_all_data(
    dist_cutoff = 50, num_primary_beta = 5, noise_level=None):
    np.random.seed(0)
    grid_size = 2000
    num_dots = 8000
    num_trivial_beta = 50 - num_primary_beta
    locations = simulate_grid(grid_size, num_dots)
    senders, receivers = simulate_senders_receivers(
        num_dots, locations,
        centers = np.array(
            [[292, 799],[894, 149],[855, 774],[611, 391],[297, 366]]
            ) * 2,
        radius_min = 150, radius_max = 300, ring_r = 200
    )
    meta_data = pd.DataFrame(locations, columns=['X','Y'])
    meta_data['Celltype'] = 'Others'
    meta_data.iloc[receivers, 2] = 'Receivers'
    meta_data.iloc[senders, 2] = 'Senders'
    pip = simulate_pip(receivers, senders, locations, dist_cutoff)
    total_ip = simulate_pip(receivers, senders, locations, dist_cutoff*1.5)

    length_checker = np.vectorize(len)
    receivers = np.array(receivers)[length_checker(total_ip)>0]
    pip = np.array(pip)[length_checker(total_ip)>0]
    total_ip = np.array(total_ip)[length_checker(total_ip)>0]

    # simulate sender expression with normaal distribution
    # This also means in the real application, the pathways/gene expression need to be
    # z-scored
    exp_sender = pd.DataFrame(
        np.random.normal(size=[len(senders), 50]), index = senders)

    # simulate betas
    primary_beta = np.random.normal(0, 1, size=num_primary_beta) * 10
    while (np.abs(primary_beta) < 10).any():
        primary_beta = [2*x if abs(x)<10 else x for x in primary_beta]
    trivial_beta = np.random.normal(0 ,1, num_trivial_beta)
    betas = np.append(primary_beta, trivial_beta)

    # simulate receiver expression
    exp_receiver = np.zeros(shape=(len(receivers), 1))
    for i in range(len(receivers)):
        i_senders = pip[i]
        exp_receiver[i] = np.matmul(exp_sender.loc[i_senders], betas).sum() + np.random.normal(-0.2, 0.1)
    exp_receiver = (exp_receiver > 0) + 0
    if noise_level is not None:
        exp_receiver = exp_receiver + np.random.normal(0, noise_level, exp_receiver.shape)
        exp_receiver = (exp_receiver > 0) + 0
    return exp_sender, exp_receiver, total_ip, pip, receivers, senders, meta_data, betas

def write_simulation_files(
    exp_sender, exp_receiver, total_ip, pip, receivers, senders, meta_data, betas, 
    log_info, output_path):
    # Write out simulation files
    # exp sender 
    exp_sender_fn = output_path + '/exp_sender.json'
    exp_sender_dict = dict()
    for i, s in enumerate(total_ip):
        exp_sender_dict['r' + str(i+1)] = exp_sender.loc[s].values.tolist()
    with open(exp_sender_fn, 'w') as fp:
        json.dump(exp_sender_dict, fp)

    # dist_sender
    dist_sender_fn = output_path + '/dist_sender.json'
    sender_dist_dict = dict()
    dist_r2s = cdist(meta_data.iloc[receivers, :2], meta_data.iloc[senders, :2])
    dist_r2s = pd.DataFrame(dist_r2s, columns = senders)
    for i, row in dist_r2s.iterrows():
        sender_dist_dict['r' + str(i+1)] = row.loc[total_ip[i]].tolist()
    with open(dist_sender_fn, 'w') as fp:
        json.dump(sender_dist_dict, fp)

    # exp receiver
    exp_receiver_fn = output_path + '/exp_receiver.csv'
    pd.DataFrame(
        exp_receiver).to_csv(exp_receiver_fn, 
        header=None, 
        index=None)

    pd.Series(betas).to_csv(output_path + '/betas.csv', header=None, index=None)

    # metadata 
    meta_data_fn = output_path + '/simulation_metadata.txt'
    meta_data['Sender_cells'] = ''
    meta_data.iloc[receivers,3] = [
        ','.join([str(i) for i in x]) for x in total_ip]
    meta_data['Sender_cells_PI'] = ''
    meta_data.iloc[receivers,4] = [
        ','.join([str(i) for i in x]) for x in pip]
    meta_data.to_csv(meta_data_fn, sep='\t')

    # write log
    with open(output_path + '/log.txt', 'w') as f:
        f.write(log_info)
    return exp_sender_fn, dist_sender_fn, exp_receiver_fn, meta_data_fn

def simulation_worker(cmd):
    # A simple wrapper for running denoise jobs.
    print(cmd)
    os.system(cmd)
    # os.system('rm -r output_fn') # Get rid of output.
    return
#%%
spacia_path = '/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/code/spacia'

try:
    dist_cutoff = sys.argv[1]
except:
    dist_cutoff = 50

try:
    num_primary_beta = sys.argv[2]
except:
    num_primary_beta = 5

try:
    noise_level = float(sys.argv[3])
except:
    noise_level = None

if noise_level is None:
    sim_prefix = 'base'
else:
    sim_prefix = 'noise_level_{:.1f}'.format(noise_level)
# setting up output folder
output_path = spacia_path.replace(
    'spacia', 'data/simulation/{}_{}_pathways'.format(sim_prefix, num_primary_beta))

#! Make sure the log info is accurate!!!
spacia_job_prefix = '{}_{}_pathways'.format(sim_prefix, num_primary_beta)
log_info = spacia_job_prefix
log_info = '\n'.join([
    log_info,
    '2000 x 2000 grid of cells, 5 blobs.',
    'True sender distance cutoff: {}.'.format(dist_cutoff),
    'Sender exp = Normal expression',
    'Number of sender pathways: {}'.format(num_primary_beta),
    'spacia job: ' + spacia_job_prefix,
    ]
)
if not os.path.exists(output_path):
    os.makedirs(output_path)
exp_sender, exp_receiver, total_ip, pip, receivers, senders, meta_data, betas = simulate_all_data(
    dist_cutoff = dist_cutoff, 
    num_primary_beta = num_primary_beta, 
    noise_level = noise_level
)
# Write out simulation files
exp_sender_fn, dist_sender_fn, exp_receiver_fn, meta_data_fn = write_simulation_files(
    exp_sender, exp_receiver, total_ip, pip, receivers, senders, meta_data, 
    betas, log_info, output_path)

# Run spacia job
spacia_R_job_path = os.path.join(spacia_path, 'spacia_job.R')
ntotal = np.arange(10000,100000,10000)
nthin = 10
nchain = 2

spacia_jobs = []
for n in ntotal:
    nwarm = int(min(20000, n/2))
    job_id = '{}_Ntotal_{}Nwarm_{}'.format(spacia_job_prefix, n, nwarm)
    spacia_output_path = os.path.join(output_path, job_id)
    if not os.path.exists(spacia_output_path):
        os.makedirs(spacia_output_path)
    spacia_jobs.append(' '.join([
        'Rscript',
        spacia_R_job_path,
        spacia_path + '/',
        exp_sender_fn,
        dist_sender_fn,
        exp_receiver_fn,
        job_id,
        str(n),
        str(nwarm),
        str(nthin),
        str(nchain),
        spacia_output_path + '/',
        ]
        )
    )

with Pool(16) as p:
    _ = p.map(simulation_worker, spacia_jobs)


'''
spacia_path = args[1]
exp_sender = args[2]
dist_receiver_sender = args[3]
receivers_mtx = args[4]
job_id = args[5]
ntotal = as.integer(args[6])
nwarm= as.integer(args[7])
nthin= as.integer(args[8])
nchain= as.integer(args[9])
output_path = args[10] # output path need to have '/' at the end
'''
# %%
# 一个是viaration of simulated data。这个你弄起来容易。
# 最好还有variation of model prior parameter distribution。
# 你要是搞得定的话，你自己试试。不行的话，我教你。就是R code里面要tweak点东西。其实也没多难