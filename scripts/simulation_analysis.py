#%%
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from collections import defaultdict
import json
import os

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

#%%
sns.set_theme('paper', font = 'Arial', font_scale=2, style='white')
simulation_input = '/project/shared/xiao_wang/projects/cell2cell_inter/code/data/simulation'
simulation_output = '/project/shared/xiao_wang/projects/cell2cell_inter/code/data/simulation/figures'
if not os.path.exists(simulation_output):
    os.makedirs(simulation_output)

#%%   
_ = plt.figure(figsize=(10,8))
sns.scatterplot(
    data = metadata, x = 'X', y = 'Y', hue='Celltype')
plt.legend(bbox_to_anchor = (1, 0.5), loc = 'center left')
plt.tight_layout()
plt.savefig(simulation_output + '/blobs.pdf')

_ = plt.figure(figsize=(8,8))
pip_n = [len(x) for x in pip]
total_ip_n = [len(x) for x in total_ip]
num_senders = pd.DataFrame(pip_n, columns = ['Num_Senders'])
num_senders['Type'] = 'Primary only'
num_senders = num_senders.append(
    pd.DataFrame(total_ip_n, columns = ['Num_Senders']))
num_senders.fillna('All Possible Senders', inplace = True)
sns.displot(
    data = num_senders, x = 'Num_Senders', hue = 'Type', multiple="stack")
plt.savefig(simulation_output + '/PIP distributions.pdf')


plot_pip(pip, receivers, metadata, figsize=(16,10))
plt.tight_layout()
plt.savefig(simulation_output + '/True interactions.pdf')
plot_pip(pip, receivers, metadata, figsize=(8,4), bbox = (500,800,200))
plt.tight_layout()
plt.savefig(simulation_output + '/True interactions zoom.pdf')

# %%
simulation_folder = output_path
job_id = ''
b = pd.read_csv(os.path.join(output_path, job_id + '_b.txt'), sep='\t')
res_beta = pd.read_csv(os.path.join(output_path, job_id + '_beta.txt'), sep='\t')
fdr = pd.read_csv(os.path.join(output_path, job_id + '_FDRs.txt'), sep='\t')
pip_res = pd.read_csv(os.path.join(output_path, job_id + '_pip.txt'), sep='\t')
# pip_res = pd.read_csv(os.path.join(output_path, job_id + '_pip_res.txt'), sep='\t')
# %%
pi_vector = []
for i, s in enumerate(total_ip):
    pi_vector += [True if x in pip[i] else False for x in total_ip[i]]
pip_res['Primary_instance'] = pi_vector
pip_res = pip_res.melt(
    value_vars = pip_res.columns[:-1], id_vars='Primary_instance', 
    value_name = 'Pi_score', var_name = 'Chain')
# %%
sns.violinplot(
    data = pip_res, x = 'Primary_instance', y = 'Pi_score', hue='Chain')
plt.legend(bbox_to_anchor = (1,0.5), loc = 'center left')
plt.tight_layout()
plt.savefig(simulation_output + '/Primary vs non-Primary pip score.pdf')

# %%
import numpy as np, scipy.stats as st
for col in res_beta.columns:
    interval = st.t.interval(
        0.95, res_beta.shape[0]-1, loc=np.mean(betas[col]), 
        scale=st.sem(betas[col]))
    print(col, '{:.3f} to {:.3f}'.format(*interval))
for b in betas:
    print(round(b,2))
# %%
true_beta = betas
all_chains_mean = res_beta.mean().to_frame(name='Correlation')
all_chains_mean['Chains'] = 'All'
for i in range(5):
    res_beta_mean = res_beta.iloc[i*6000:(i+1)*6000].mean()
    all_chains_mean = all_chains_mean.append(
        res_beta_mean.to_frame(name = 'Correlation'))
    all_chains_mean.fillna('Chain' + str(i+1), inplace=True)
all_chains_mean['True beta'] = true_beta.tolist()*6
cor = np.corrcoef(true_beta, all_chains_mean.iloc[:50,0])[0,1]
_ = plt.figure(figsize=(8,6))
sns.scatterplot(
    data = all_chains_mean, x = 'Correlation', y = 'True beta', hue = 'Chains',
    s = 25, alpha = 0.5, linewidth = 0)
plt.legend(bbox_to_anchor = (1,0.5), loc = 'center left')
plt.title('Correlation between True and Spacia Beta across all chains: {:.3f}'.format(cor))
plt.tight_layout()
plt.savefig(simulation_output + '/Beta per chain correlation.pdf')
# %%
