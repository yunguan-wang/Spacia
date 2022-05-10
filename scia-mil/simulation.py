#%%
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from collections import defaultdict
import json
import os

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

def plot_pip(pip, receivers, locations):
    pip_all = np.unique([x for p in pip for x in p])
    receivers_all = [receivers[i] for i in range(len(receivers)) if pip[i].shape[0]>0]
    _ = plt.figure(figsize=(8,8))
    sns.scatterplot(
        x = locations[:,0], y = locations[:,1], color='gray', s=10,
        label = 'Others')
    sns.scatterplot(
        x = locations[receivers,0], y = locations[receivers,1], 
        color = 'g', s = 10, label = 'Out-of-range Receivers')
    sns.scatterplot(
        x = locations[receivers_all,0], y = locations[receivers_all,1], 
        color = 'g', s = 25, label = 'Receivers')
    sns.scatterplot(
        x = locations[senders,0], y = locations[senders,1], 
        color = 'r', s = 10, label = 'Non-PIP')
    sns.scatterplot(
        x = locations[pip_all,0], y = locations[pip_all,1], color = 'r', s=50,
        label = 'PIP')
    plt.legend(bbox_to_anchor = (1,0.5), loc = 'center left')
    for i in range(len(receivers)):
        if len(pip[i]) == 0:
            continue
        else:
            r = receivers[i]
            p = pip[i]
            for indiv_p in p:
                start = locations[indiv_p]
                end = locations[r]
                dx, dy = end-start
                plt.arrow(
                    *start, dx, dy,
                    head_width = 10,
                    length_includes_head=True)

def simulate_senders_receivers(
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

np.random.seed(0)
output_path = '/project/shared/xiao_wang/projects/cell2cell_inter/simulation'
grid_size = 2000
num_dots = 8000
dist_cutoff = 50
locations = simulate_grid(grid_size, num_dots)
senders, receivers = simulate_senders_receivers(
    centers = np.array(
        [[292, 799],[894, 149],[855, 774],[611, 391],[297, 366]]
        ) * 2,
    radius_min = 150, radius_max = 300, ring_r = 100
)
pip = simulate_pip(receivers, senders, locations, dist_cutoff)
total_ip = simulate_pip(receivers, senders, locations, dist_cutoff*1.5)

length_checker = np.vectorize(len)
receivers = np.array(receivers)[length_checker(total_ip)>0]
pip = np.array(pip)[length_checker(total_ip)>0]
total_ip = np.array(total_ip)[length_checker(total_ip)>0]

# simulate sender expression with normaal distribution
# This also means in the real application, the pathways/gene expression need to be
# z-scored
sender_names = ['s' + str(x) for x in senders]
exp_sender = pd.DataFrame(
    np.random.normal(size=[len(senders), 50]), index = senders)

# simulate betas
primary_beta = np.random.normal(0, 1, size=10) * 10
while (np.abs(primary_beta) < 10).any():
    primary_beta = [2*x if abs(x)<10 else x for x in primary_beta]
trivial_beta = np.random.normal(0 ,1, 40)
betas = np.append(primary_beta, trivial_beta)

# simulate receiver expression
exp_receiver = np.zeros(shape=(len(receivers), 1))
first = []
second = []
for i in range(len(receivers)):
    i_senders = pip[i]
    first.append(np.matmul(exp_sender.loc[i_senders], betas).sum())
    second.append(np.random.normal(-0.2, 0.1))
    exp_receiver[i] = np.matmul(exp_sender.loc[i_senders], betas).sum() + np.random.normal(-0.2, 0.1)
exp_receiver = (exp_receiver > 0) + 0

# # construct distance matrix
# dist_r2s = cdist(locations[receivers], locations[senders])
# dist_r2s = pd.DataFrame(dist_r2s, columns = sender_names)
# dist_r2s.to_csv(output_path + '/dist_r2s.csv')

# # construct sender exp matrix
# exp_sender.index = sender_names
# exp_sender.to_csv(output_path + '/exp_sender.csv')

# # construct receiver exp matrix
# receiver_df = pd.DataFrame(total_ip.copy(), columns = ['senders'])
# receiver_df.senders = receiver_df.senders.apply(
#     lambda x: ','.join(['s' + str(i) for i in x])
#     )
# receiver_df['exp_receiver'] = exp_receiver
# receiver_df.to_csv(output_path + '/exp_receiver.csv')
# pd.Series(betas).to_csv(output_path + '/betas.csv', header=None, index=None)

exp_sender_dict = dict()
for i, s in enumerate(total_ip):
    exp_sender_dict['r' + str(i+1)] = exp_sender.loc[s].values.tolist()
with open(output_path + '/exp_sender.json', 'w') as fp:
    json.dump(exp_sender_dict, fp)

sender_dist_dict = dict()
dist_r2s = cdist(locations[receivers], locations[senders])
dist_r2s = pd.DataFrame(dist_r2s, columns = senders)
for i, row in dist_r2s.iterrows():
    sender_dist_dict['r' + str(i+1)] = row.loc[total_ip[i]].tolist()
with open(output_path + '/exp_dist.json', 'w') as fp:
    json.dump(sender_dist_dict, fp)

pd.DataFrame(
    exp_receiver).to_csv(output_path + '/exp_receiver.csv', 
    header=None, 
    index=None)

pd.Series(betas).to_csv(output_path + '/betas.csv', header=None, index=None)
#%%
# Run 
#%%
# _ = plt.figure(figsize=(8,8))
# sns.scatterplot(x = locations[:,0], y = locations[:,1], color='gray')
# sns.scatterplot(x = locations[senders,0], y = locations[senders,1], color = 'r')
# sns.scatterplot(x = locations[receivers,0], y = locations[receivers,1], color = 'g')

# _ = plt.figure(figsize=(8,8))
# pip_n = [len(x) for x in pip]
# total_ip_n = [len(x) for x in total_ip]
# plt.hist(pip_n, label = 'Pip')
# plt.hist(total_ip_n, label = 'Total_ip', alpha=0.5)
# plt.legend()
# plot_pip(pip, receivers, locations)

# _ = plt.figure(figsize=(8,8))
# plt.hist(first, label = 'primary beta')
# plt.hist(second, label = 'secondary beta')
# plt.legend()
# %%
simulation_folder = output_path
job_id = ''
b = pd.read_csv(os.path.join(output_path, job_id + '_b.txt'), sep='\t')
beta = pd.read_csv(os.path.join(output_path, job_id + '_beta.txt'), sep='\t')
fdr = pd.read_csv(os.path.join(output_path, job_id + '_FDRs.txt'), sep='\t')
pip_res = pd.read_csv(os.path.join(output_path, job_id + '_pip.txt'), sep='\t')
# pip_res = pd.read_csv(os.path.join(output_path, job_id + '_pip_res.txt'), sep='\t')
# %%
pi_vector = []
for i, s in enumerate(total_ip):
    pi_vector += [True if x in pip[i] else False for x in total_ip[i]]
# %%
pip_res['Primary_instance'] = pi_vector
pip_res = pip_res.melt(
    value_vars = pip_res.columns[:4], id_vars='Primary_instance', 
    value_name = 'Pi_score', var_name = 'Chain')
# %%
sns.violinplot(data = pip_res, x = 'Primary_instance', y = 'Pi_score', hue='Chain')
# %%
sns.dis(data = pip_res, x = 'Primary_instance', y = 'pip_recal')
# %%
