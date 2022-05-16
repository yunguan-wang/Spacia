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

def plot_pip(
    s_list, r_list, metadata, figsize = (10,8),
    bbox = None # (500, 800, 200) as X, Y for top left corner and side
    ):
    pip_all = np.unique([x for p in s_list for x in p])
    receivers_all = [r_list[i] for i in range(len(r_list)) if s_list[i].shape[0]>0]
    plot_data = metadata.copy()
    plot_data['dot_size'] = 10
    plot_data.Celltype.replace(
        'Receivers','Out-of-range Receiver Cell Type', inplace = True)
    plot_data.Celltype.replace(
        'Senders','Out-of-range Sender Cell Type', inplace = True)
    plot_data.iloc[pip_all,2] = 'Sender Cell'
    plot_data.iloc[receivers_all,2] = 'Receiver Cell'
    plot_data.iloc[pip_all,5] = 20
    plot_data.iloc[receivers_all,5] = 20
    if bbox is not None:
        X_min, Y_min, delta = bbox
        mask = (plot_data.X >= X_min) & (plot_data.X <= X_min + delta)
        mask = mask & (plot_data.Y >= Y_min) & (plot_data.Y <= Y_min + delta)
        plot_data = plot_data[mask]
        receivers_mask = [
            True if x in plot_data.index else False for x in r_list]
        r_list = np.array(r_list)[receivers_mask]
        s_list = np.array(s_list)[receivers_mask]
        plot_data.dot_size = plot_data.dot_size * 10

    _ = plt.figure(figsize=figsize)
    sns.scatterplot(
        data = plot_data, x = 'X', y = 'Y', hue='Celltype', linewidth=0,
        size = 'dot_size',
        hue_order = [
            'Others', 'Out-of-range Receiver Cell Type', 'Receiver Cell',
       'Out-of-range Sender Cell Type', 'Sender Cell'],
       palette = ['grey', 'g', 'g', 'r', 'r']
       )
    plt.legend(bbox_to_anchor = (1,0.5), loc = 'center left')

    for i in range(len(r_list)):
        if len(s_list[i]) == 0:
            continue
        else:
            r = r_list[i]
            p = s_list[i]
            for indiv_p in p:
                try:
                    start = plot_data.loc[indiv_p,:'Y']
                    end = plot_data.loc[r,:'Y']
                except KeyError:
                    continue
                dx, dy = end-start
                plt.arrow(
                    *start, dx, dy,
                    head_width = 5,
                    length_includes_head=True, 
                    edgecolor=None,
                    color = 'k')

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

def write_simulation_files(
    exp_sender, total_ip, receivers, senders, metadata, betas, output_path):
    # Write out simulation files
    exp_sender_dict = dict()
    for i, s in enumerate(total_ip):
        exp_sender_dict['r' + str(i+1)] = exp_sender.loc[s].values.tolist()
    with open(output_path + '/exp_sender.json', 'w') as fp:
        json.dump(exp_sender_dict, fp)

    sender_dist_dict = dict()
    dist_r2s = cdist(metadata.iloc[receivers, :2], metadata.iloc[senders, :2])
    dist_r2s = pd.DataFrame(dist_r2s, columns = senders)
    for i, row in dist_r2s.iterrows():
        sender_dist_dict['r' + str(i+1)] = row.loc[total_ip[i]].tolist()
    with open(output_path + '/dist_sender.json', 'w') as fp:
        json.dump(sender_dist_dict, fp)

    pd.DataFrame(
        exp_receiver).to_csv(output_path + '/exp_receiver.csv', 
        header=None, 
        index=None)

    pd.Series(betas).to_csv(output_path + '/betas.csv', header=None, index=None)

    metadata['Sender_cells'] = ''
    metadata.iloc[receivers,3] = [
        ','.join([str(i) for i in x]) for x in total_ip]
    metadata['Sender_cells_PI'] = ''
    metadata.iloc[receivers,4] = [
        ','.join([str(i) for i in x]) for x in pip]
    metadata.to_csv(output_path + '/simulation_metadata.txt', sep='\t')
#%%
np.random.seed(0)
spacia_path = '/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/code/spacia'
output_path = spacia_path.replace('spacia', 'data/simulation/base')
if not os.path.exists(output_path):
    os.makedirs(output_path)
grid_size = 2000
num_dots = 8000
dist_cutoff = 50
locations = simulate_grid(grid_size, num_dots)
senders, receivers = simulate_senders_receivers(
    centers = np.array(
        [[292, 799],[894, 149],[855, 774],[611, 391],[297, 366]]
        ) * 2,
    radius_min = 150, radius_max = 300, ring_r = 200
)
metadata = pd.DataFrame(locations, columns=['X','Y'])
metadata['Celltype'] = 'Others'
metadata.iloc[receivers, 2] = 'Receivers'
metadata.iloc[senders, 2] = 'Senders'
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
for i in range(len(receivers)):
    i_senders = pip[i]
    exp_receiver[i] = np.matmul(exp_sender.loc[i_senders], betas).sum() + np.random.normal(-0.2, 0.1)
exp_receiver = (exp_receiver > 0) + 0

# Write out simulation files
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
with open(output_path + '/dist_sender.json', 'w') as fp:
    json.dump(sender_dist_dict, fp)

pd.DataFrame(
    exp_receiver).to_csv(output_path + '/exp_receiver.csv', 
    header=None, 
    index=None)

pd.Series(betas).to_csv(output_path + '/betas.csv', header=None, index=None)

metadata['Sender_cells'] = ''
metadata.iloc[receivers,3] = [
    ','.join([str(i) for i in x]) for x in total_ip]
metadata['Sender_cells_PI'] = ''
metadata.iloc[receivers,4] = [
    ','.join([str(i) for i in x]) for x in pip]
metadata.to_csv(output_path + '/simulation_metadata.txt', sep='\t')


# %%
# 一个是viaration of simulated data。这个你弄起来容易。
# 最好还有variation of model prior parameter distribution。你要是搞得定的话，你自己试试。不行的话，我教你。就是R code里面要tweak点东西。其实也没多难