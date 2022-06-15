#%%
import pandas as pd
from requests import head
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from collections import defaultdict
import json
import os
from sklearn.metrics import roc_curve, auc

def plot_pip(
    s_list, r_list, meta_data, figsize = (10,8),
    bbox = None, # (500, 800, 200) as X, Y for top left corner and side
    **kwargs
    ):
    pip_all = np.unique([x for p in s_list for x in p])
    receivers_all = [r_list[i] for i in range(len(r_list)) if len(s_list[i])>0]
    plot_data = meta_data.copy()
    plot_data['dot_size'] = 10
    plot_data.Celltype.replace(
        'Receivers','Out-of-range Receiver Cell Type', inplace = True)
    plot_data.Celltype.replace(
        'Senders','Out-of-range Sender Cell Type', inplace = True)
    plot_data.iloc[pip_all,2] = 'Sender Cell'
    plot_data.iloc[receivers_all,2] = 'Receiver Cell'
    plot_data.iloc[pip_all,5] = 50
    plot_data.iloc[receivers_all,5] = 50
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
        sizes = (25,150),
        hue_order = [
            'Others', 'Out-of-range Receiver Cell Type', 'Receiver Cell',
       'Out-of-range Sender Cell Type', 'Sender Cell'],
       **kwargs
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
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["legend.markerscale"] = 2
# simulation_input = '/project/shared/xiao_wang/projects/cell2cell_inter/data/nthin100/base_10_pathways/spacia_base_10_pathways_Ntotal_80000_Nwarm_20000_Nthin_100'
simulation_input = '/project/shared/xiao_wang/projects/cell2cell_inter/data/simulation/base_5chains_nthin10_ntotal80000'
simulation_output = '/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/simulation_figures'
if not os.path.exists(simulation_output):
    os.makedirs(simulation_output)
#%%
# Choose colors
colors = sns.color_palette("bright").as_hex()
#%%
# ==============================================================================
# simulation data specs
meta_data = pd.read_csv(
    os.path.join(simulation_input,'simulation_metadata.txt'), index_col=0, sep='\t')
receivers = meta_data.dropna(subset = ['Sender_cells']).index
candidate_senders = meta_data.loc[receivers,'Sender_cells'].str.split(',').tolist()
true_senders = meta_data.loc[receivers,'Sender_cells_PI'].str.split(',').tolist()
true_senders = [[] if isinstance(x, float) else x for x in true_senders]
candidate_senders = [list(map(lambda x: int(x), y)) for y in candidate_senders]
true_senders = [list(map(lambda x: int(x), y)) for y in true_senders]
# job_id = 'spacia_base_5_pathways_Ntotal_50000_Nwarm_20000_Nthin_100'
job_id = ''
pred_b = pd.read_csv(os.path.join(simulation_input, job_id + '_b.txt'), sep='\t')
# true_beta = pd.read_csv(os.path.join(simulation_input, 'betas.csv'))
true_beta = pd.read_csv(os.path.join(simulation_input, 'betas.csv'), header=None)
res_beta = pd.read_csv(os.path.join(simulation_input, job_id + '_beta.txt'), sep='\t')
fdr = pd.read_csv(os.path.join(simulation_input, job_id + '_FDRs.txt'), sep='\t')
pip_res = pd.read_csv(os.path.join(simulation_input, job_id + '_pip.txt'), sep='\t')

# %%
# ==============================================================================
# spacia results
pi_vector = []
for i, s in enumerate(candidate_senders):
    pi_vector += [True if x in true_senders[i] else False for x in candidate_senders[i]]
pip_res['Primary_instance'] = pi_vector
plotdata = pip_res.melt(
    value_vars = pip_res.columns[:-1], id_vars='Primary_instance', 
    value_name = 'Pi_score', var_name = 'Chain')
plotdata.Chain = plotdata.Chain.str.replace('pip.', 'Chain ')
_ = plt.figure(figsize=(8,6))
sns.violinplot(
    data = plotdata, x = 'Primary_instance', y = 'Pi_score', hue='Chain',
    palette = colors[:5])
plt.xlabel('')
plt.ylabel('Primary instance score')
plt.legend(bbox_to_anchor = (1,0.5), loc = 'center left')
plt.tight_layout()
plt.savefig(simulation_output + '/Primary vs non-Primary true_senders score.pdf')
#%%
fpr = dict()
tpr = dict()
roc_auc = dict()
for i, chain in enumerate(plotdata.Chain.unique()):
    pip_chain = plotdata[plotdata.Chain==chain]
    fpr[i], tpr[i], _ = roc_curve(pip_chain.Primary_instance, pip_chain.Pi_score)
    roc_auc[i] = auc(fpr[i], tpr[i])

plt.figure(figsize=(8,4))
lw = 5
for i, chain in enumerate(plotdata.Chain.unique()):
    plt.plot(
        fpr[i],
        tpr[i],
        color=colors[i],
        lw=lw,
        label= "{} (area = {:.2f})".format(chain, roc_auc[i]),
    )
plt.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend(bbox_to_anchor = (1, 0.5), loc="center left")
plt.tight_layout()
plt.savefig(os.path.join(simulation_output, 'True senders ROC.pdf'))
plt.close()
# %%
import numpy as np, scipy.stats as st
for col in res_beta.columns:
    interval = st.t.interval(
        0.95, res_beta.shape[0]-1, loc=np.mean(res_beta[col]), 
        scale=st.sem(res_beta[col]))
    print(col, '{:.3f} to {:.3f}'.format(*interval))
# %%
all_chains_mean = res_beta.mean().to_frame(name='Predicted')
all_chains_mean['Chains'] = 'All'
for i in range(5):
    res_beta_mean = res_beta.iloc[i*6000:(i+1)*6000].mean()
    all_chains_mean = all_chains_mean.append(
        res_beta_mean.to_frame(name = 'Predicted'))
    all_chains_mean.fillna('Chain' + str(i+1), inplace=True)
all_chains_mean['True beta'] = true_beta.iloc[:,0].tolist()*6
cor = np.corrcoef(true_beta.values.reshape(50,), all_chains_mean.iloc[:50,0])[0,1]
_ = plt.figure(figsize=(8,6))
# sns.scatterplot(
#     data = all_chains_mean, x = 'Predicted', y = 'True beta', hue = 'Chains',
#     s = 50, alpha = 0.8, linewidth = 0)
plotdata = all_chains_mean[all_chains_mean.Chains == 'All']
plotdata['Beta type'] = ['True'] * 10 + ['Trivial'] * 40
sns.scatterplot(
    data = plotdata, x = 'Predicted', y = 'True beta', hue = 'Beta type',
    s = 100, alpha = 0.8, linewidth = 1, palette='bright')
plt.legend(
    bbox_to_anchor = (1,0.5), loc = 'center left', markerscale=2)
plt.title('Correlation between True and Spacia Beta across all chains: {:.3f}'.format(cor))
plt.xlabel('Predicted Beta')
plt.tight_layout()
plt.savefig(simulation_output + '/Beta correlation.pdf')
# %%
beta_ranks = pd.DataFrame(
    np.arange(1,51,1),
    index = [
        'beta.' + str(i+1) for i in (abs(true_beta.iloc[:,0]).sort_values(ascending=False)).index
        ],
    columns = ['True'])
predicted_ranks = res_beta.mean().to_frame(name = 'Predicted')
beta_ranks['Predicted'] = np.argsort(abs(predicted_ranks.Predicted.values))[::-1] + 1
#%%
non_trivial_betas = ['beta.' + str(i+1) for i in range(5)]
true_beta_ranks = beta_ranks[beta_ranks.index.isin(non_trivial_betas)]
trivial_beta_ranks = beta_ranks[~beta_ranks.index.isin(non_trivial_betas)]
_ = plt.figure(figsize=(4,7))
for i, ranks in enumerate([true_beta_ranks, trivial_beta_ranks]):
    nrows = ranks.shape[0]
    x = [0]*nrows + [1]*nrows
    y = ranks.values.T.flatten().tolist()
    plt.scatter(
        x, y, c = ['red', 'grey'][i], 
        label = ['Significant Pathway', 'Trivial Pathway'][i])
for _, row in true_beta_ranks.iterrows():
    x = 0
    y = row[0]
    dx = 1
    dy = row[1]-row[0]
    plt.arrow(
        x, y, dx, dy,
        head_width = 0,
        length_includes_head=True, 
        edgecolor=None,
        color = 'k')
plt.xticks([0,1],['True Rank', 'Predicted Rank'])
plt.yticks([],[])
plt.ylabel('Beta Ranks of Sender Pathways')
plt.legend(bbox_to_anchor = (0.5, -0.2), loc = 'center')
plt.tight_layout()
plt.savefig(simulation_output + '/Beta ranks.pdf')
#%%   
_ = plt.figure(figsize=(10,8))
sns.scatterplot(
    data = meta_data, x = 'X', y = 'Y', hue='Celltype', palette=colors[:3])
plt.legend(bbox_to_anchor = (1, 0.5), loc = 'center left')
plt.tight_layout()
plt.savefig(simulation_output + '/blobs.pdf')

_ = plt.figure(figsize=(8,8))
pip_n = [len(x) for x in true_senders]
c_senders_n = [len(x) for x in candidate_senders]
num_senders = pd.DataFrame(pip_n, columns = ['Num_Senders'])
num_senders['Type'] = 'Primary only'
num_senders = num_senders.append(
    pd.DataFrame(c_senders_n, columns = ['Num_Senders']))
num_senders.fillna('All Possible Senders', inplace = True)
sns.displot(
    data = num_senders, x = 'Num_Senders', hue = 'Type', multiple="stack",
    palette = 'bright')
plt.savefig(simulation_output + '/PIP distributions.pdf')


plot_colors = [colors[i] for i in [7,9,0,6,3]]
plot_pip(
    true_senders, receivers, meta_data, figsize=(16,10),
    palette=plot_colors)
plt.tight_layout()
plt.savefig(simulation_output + '/True interactions.pdf')
plot_pip(
    true_senders, receivers, meta_data, figsize=(8,4), bbox = (500,800,200),
    palette=plot_colors)
plt.tight_layout()
plt.savefig(simulation_output + '/True interactions zoom.pdf')
# %%
## Simulation plot evaluating effect of nTotal
simulation_input = '/project/shared/xiao_wang/projects/cell2cell_inter/data/nthin100/base_10_pathways'
simulation_output = '/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/simulation_figures'
if not os.path.exists(simulation_output):
    os.makedirs(simulation_output)
folders = [x for x in os.listdir(simulation_input) if os.path.isdir(simulation_input + '/' +x)]
fpr = dict()
tpr = dict()
roc_auc = dict()
labels = []
for k, fd in enumerate(sorted(folders)):
    labels.append(fd.split('_')[5])
    plotdata = pd.read_csv(
        os.path.join(simulation_input,fd, fd + '_pip.txt'), sep='\t').mean(axis=1).to_frame(name = 'PI_score')
    pi_vector = []
    for i, s in enumerate(candidate_senders):
        pi_vector += [True if x in true_senders[i] else False for x in candidate_senders[i]]
    plotdata['Primary_instance'] = pi_vector

    fpr[k], tpr[k], _ = roc_curve(plotdata.Primary_instance, plotdata.PI_score)
    roc_auc[k] = auc(fpr[k], tpr[k])


plt.figure(figsize=(8,4))
lw = 2
for i, label in enumerate(labels):
    plt.plot(
        fpr[i],
        tpr[i],
        color=colors[i],
        lw=lw,
        label= 'Total MCMC chains:' + str(label),
    )
plt.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
plt.xlim([-0.01, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend(bbox_to_anchor = (1, 0.5), loc="center left")
plt.tight_layout()
plt.savefig(os.path.join(simulation_output, 'True senders ROC ntotal performance.pdf'))
plt.close()
# %%
simulation_input = '/project/shared/xiao_wang/projects/cell2cell_inter/data/nthin100'
simulation_output = '/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/simulation_figures'
if not os.path.exists(simulation_output):
    os.makedirs(simulation_output)
folders = [
    x for x in os.listdir(simulation_input) if (os.path.isdir(simulation_input + '/' +x)) &
    ('noise' in x)]
fpr = dict()
tpr = dict()
roc_auc = dict()
labels = []
for k, fd in enumerate(sorted(folders)):
    subfd = 'spacia_' + fd + '_Ntotal_50000_Nwarm_20000_Nthin_100'
    labels.append(fd.split('_')[2])
    plotdata = pd.read_csv(
        os.path.join(
            simulation_input,fd, subfd, subfd + '_pip.txt'
            ), sep='\t').mean(axis=1).to_frame(name = 'PI_score')
    pi_vector = []
    for i, s in enumerate(candidate_senders):
        pi_vector += [True if x in true_senders[i] else False for x in candidate_senders[i]]
    plotdata['Primary_instance'] = pi_vector

    fpr[k], tpr[k], _ = roc_curve(plotdata.Primary_instance, plotdata.PI_score)
    roc_auc[k] = auc(fpr[k], tpr[k])


plt.figure(figsize=(8,4))
lw = 2
for i, label in enumerate(labels):
    plt.plot(
        fpr[i],
        tpr[i],
        color=colors[i],
        lw=lw,
        label= 'Noise level:' + str(label),
    )
plt.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
plt.xlim([-0.01, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend(bbox_to_anchor = (1, 0.5), loc="center left")
plt.tight_layout()
plt.savefig(os.path.join(simulation_output, 'True senders ROC noise performance.pdf'))
plt.close()# %%

# %%
