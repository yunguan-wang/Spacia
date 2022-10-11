#%%
import pandas as pd
import scanpy as sc
import numpy as np
import os
from sklearn.mixture import GaussianMixture
from sklearn.svm import SVC
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from scipy.spatial.distance import cdist
import seaborn as sns
import matplotlib.ticker as ticker
from skimage import morphology, io
#%%
# Cycif
input_path = 'C:/Users/Jake/Downloads/'
output_path = 'C:/Users/Jake/Downloads/'
sc.settings.figdir = output_path
sns.set_theme(
    context='paper', style='white', palette='gist_rainbow', 
    font= 'Arial',font_scale=3)
sns.set_style({'ytick.left': True, 'xtick.bottom': True})
# %%
cpm = pd.read_csv(input_path + 'LUNG-1-LN_40X_master.csv')
cpm.index = ['cell_' + str(x) for x in cpm.index]
cpm = cpm[np.isfinite(cpm).all(axis=1)]
cpm_exp = cpm[
    [
        'DAPI9', 'LAG3', 'ARL13B', 'DAPI4', 'KI67', 'KERATIN', 'PD1', 'CD45RB','CD3D', 'PDL1',
        'CD4', 'CD45', 'CD8A', 'CD163','CD68', 'CD14', 'CD11B', 'FOXP3', 'CD21', 'IBA1','ASMA',
        'CD20', 'CD19', 'GFAP', 'GTUBULIN','LAMINAC', 'BANF1', 'LAMINB']]
cpm_meta = cpm[
    [x for x in cpm.columns if x not in cpm_exp.columns]
].copy()
for c in ['A488', 'A555','A647']:
    b_col = [x for x in cpm_meta.columns if c + 'background' in x]
    cpm_meta[c] = cpm_meta[b_col].mean(axis=1).values

def norm_cycif(cpm, cpm_meta, c1, c0):
    return cpm[c1] - cpm_meta.loc[cpm.index, c0]
# %%
fov = [1,1]
cells = cpm_meta.index[(cpm_meta.iloc[:,0] == fov[0]) & (cpm_meta.iloc[:,1] == fov[1])]
fov_cpm = cpm_exp.loc[cells]
fov_meta = cpm_meta.loc[cells]
# y = fov_meta['Y_position'].copy()
fov_meta['Y_position'] = fov_meta['Y_position'].max() - fov_meta['Y_position']
# fov_meta['X_position'] = y
fov_cpm = sc.AnnData(fov_cpm)
sc.tl.pca(fov_cpm)
sc.pp.neighbors(fov_cpm)
sc.tl.leiden(fov_cpm)
sc.tl.umap(fov_cpm)
fov_cpm = fov_cpm[fov_cpm.obs.leiden.isin([
    '0', '2', '9', '11', '12', '13', '14', '15', '16'])]
# %%
sc.pl.umap(fov_cpm, color = ['leiden'], legend_loc = 'on data')
sc.pl.umap(fov_cpm, color = fov_cpm.var_names, ncols=2, cmap='viridis',vmin='p5',vmax='p95')

# %%
# B/Mono/DC
# c_s = fov_cpm.obs_names[fov_cpm.obs.leiden.isin(['14','15','16','2'])].tolist()
# c_r = fov_cpm.obs_names[fov_cpm.obs.leiden.isin(['0','13'])].tolist()
# m2 = fov_cpm.obs_names[
#     (fov_cpm.obs.leiden=='2') & 
#     (fov_cpm.to_df()['CD68'] > 6.5) & 
#     (fov_cpm.to_df()['CD163'] > 9)].tolist()
# c_s = [x for x in c_s if x not in m2]
# B/T
c_s = fov_cpm.obs_names[fov_cpm.obs.leiden.isin(['13'])].tolist()
c_r = fov_cpm.obs_names[fov_cpm.obs.leiden.isin(['14'])].tolist()
# %%
gene = 'CD11B'
_ = plt.figure(figsize=(10,10))
cmap = 'viridis'
ax = sns.scatterplot(
    data = fov_meta.loc[c_s+c_r], x = 'X_position', y = 'Y_position',
    hue = fov_cpm.obs.loc[c_s+c_r, 'leiden'] == '2',
    linewidth = 0, s=20, 
    # palette = cmap,
    # hue = fov_cpm.to_df().loc[c_s+c_r, gene].values, 
    # legend=None,
)
plt.legend(
    # labels=['CD8', 'M2'],
    bbox_to_anchor = (1,.5), loc='center left', markerscale=4)
plt.tight_layout()
plt.xlabel('')
plt.ylabel('')
# norm = plt.Normalize(
#     fov_cpm.to_df().loc[c_s+c_r,gene].min(), fov_cpm.to_df().loc[c_s+c_r,gene].max())
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])
# ax.figure.colorbar(sm)
# %%
# s_gene = 'CD68'
max_d = 100
xcol, ycol = ['X_position','Y_position']
cpm = fov_cpm.to_df().copy()
sender = fov_meta.loc[c_s]
receiver = fov_meta.loc[c_r]
s2r_dist = cdist(
    sender[[xcol, ycol]], receiver[[xcol, ycol]])
sender_cells = sender.index[(s2r_dist<=max_d).any(axis=1)]
receiver_cells = receiver.index
coords = fov_meta[[xcol, ycol]]

dists = cdist(coords.loc[sender_cells], coords.loc[receiver_cells])
crit1 = (dists <= max_d).any(axis=0)
crit2 = (dists <= max_d).any(axis=1)
dists = dists[crit2,:][:,crit1]
dists = pd.DataFrame(
    dists, index = sender_cells[crit2], columns = receiver_cells[crit1])
dists = dists.stack().reset_index()
dists.columns = ['sender', 'receiver', 'distance']
dists = dists[dists.distance<=max_d]
dists['Distance_bin'] = pd.cut(
    dists.distance, bins=[0, 30, max_d], include_lowest = True)
dists = dists.dropna()

s_gene = 'PD1'
r_gene = 'CD20'
plot_data = dists.copy()
# plot_data[s_gene] = norm_cycif(cpm.loc[plot_data.sender.values], fov_meta, s_gene, 'A488').values
# plot_data[r_gene] = norm_cycif(cpm.loc[plot_data.receiver.values], fov_meta, r_gene, 'A488').values
plot_data[s_gene] = cpm.loc[plot_data.sender.values,s_gene].values
plot_data[r_gene] = cpm.loc[plot_data.receiver.values,r_gene].values
# gene1sum = plot_data.groupby(
#     ['Distance_bin', 'receiver']
#     )[s_gene].mean().reset_index()
gene1sum = plot_data.groupby(
    ['Distance_bin', 'receiver']
    )[s_gene].count().reset_index()
gene1sum[s_gene] = gene1sum[s_gene].fillna(0)
sender_gene_col = s_gene + '_sum'
gene1sum.columns = ['Distance_bin', 'receiver'] + [sender_gene_col]
plot_data = plot_data.merge(gene1sum, on=['Distance_bin', 'receiver'])
plot_data = plot_data.drop_duplicates(subset = ['Distance_bin','receiver'])
# plot_data[sender_gene_col] = [
#     'High' if x>=10 else 'Low' for x in plot_data[sender_gene_col].values]
for _ ,df in plot_data.groupby('Distance_bin'):
    cuts = pd.cut(
        df[s_gene + '_sum'],bins=2,
        labels = ['Low','High'])
    plot_data.loc[
        cuts.index, sender_gene_col
        ] = cuts.values.astype(str)

cmap = ['r','cyan']
_ = plt.figure(figsize=(10,6))
ax = sns.boxplot(
    data = plot_data,
    x = 'Distance_bin',
    y = r_gene,
    hue = sender_gene_col,
    hue_order=['Low','High'],
    palette = cmap,
    linewidth=5)
plt.legend(
    bbox_to_anchor = (1, 0.5), loc = 'center left',
    title= 'Number of CD4 T cells')
plt.ylabel('Receiver ' + r_gene + ' expression', fontdict={'weight':'bold'})
plt.xlabel('')
plt.xticks(
    [0,1], labels = ['Adjacent', 'Distant'], fontweight='bold')
plt.tight_layout()
plt.savefig(output_path + '/CD4_promote_Bcells_{}.pdf'.format(r_gene))
# plt.close()
# %%
_ = plt.figure(figsize=(12,6))
pimg_c_s = plot_data.sender.unique().tolist()
pimg_c_r = plot_data.receiver.unique().tolist()
# pimg_c_s = plot_data[plot_data.distance<=30].sender.unique().tolist()
# pimg_c_r = plot_data[plot_data.distance<=30].receiver.unique().tolist()
pseudo_img = fov_meta.loc[pimg_c_s + pimg_c_r].copy()
pseudo_img['Celltype'] = fov_cpm.obs.loc[pimg_c_s + pimg_c_r, 'leiden'].astype(str).values
pseudo_img['Celltype'] = ['B cells' if x == '14' else 'CD4_Th' for x in pseudo_img['Celltype']]
pseudo_img['KI67'] = cpm.loc[pimg_c_s + pimg_c_r, 'KI67'].values
pseudo_img.loc[pimg_c_s, 'KI67'] = cpm.loc[pimg_c_s + pimg_c_r, 'KI67'].min()

# pseudo_img = pseudo_img[
#     (pseudo_img.X_position>=4000) &
#     (pseudo_img.X_position<=6000) &
#     (pseudo_img.Y_position>=3500) &
#     (pseudo_img.Y_position<=5500)]
_ = plt.figure(figsize = (12,8))
sns.scatterplot(
    data = pseudo_img[pseudo_img['Celltype']!='B cells'], 
    x = 'X_position', y = 'Y_position', linewidth = 0, s=50, alpha=0.75,
    c = ['cyan']
)
sns.scatterplot(
    data = pseudo_img[pseudo_img['Celltype']=='B cells'], 
    x = 'X_position', y = 'Y_position',
    hue = 'KI67', linewidth = 0, s=50, palette = 'Reds', alpha=0.75,
)
plt.legend(bbox_to_anchor = (1,.5), loc='center left', markerscale=4, title='KI67')
interactions = plot_data[plot_data.distance<=30]
for _, row in interactions.iterrows():
    s, r = row[:2]
    if (s in pseudo_img.index) & (r in pseudo_img.index):
        x, y = pseudo_img.loc[s, ['X_position', 'Y_position']]
        x1, y1 = pseudo_img.loc[r, ['X_position', 'Y_position']]
        dx, dy = x1-x, y1-y
        plt.arrow(x, y, dx, dy, head_width=0, head_length=0, color = 'k', lw=2)
plt.tight_layout()
plt.xlabel('')
plt.ylabel('')
plt.savefig(output_path + 'T_B_interactions.pdf')

#%%
from sklearn.preprocessing import minmax_scale
_ = plt.figure(figsize = (12,8))
pseudo_img['KI67'] = minmax_scale(pseudo_img['KI67']) + 0.01
ax = sns.kdeplot(
    data=pseudo_img, x="X_position", y="Y_position", hue="Celltype", fill=True,
    weights = 'KI67', 
)

# %%
gene = 'CD8A'
_ = plt.figure(figsize=(10,10))
cmap = 'viridis'
ax = sns.scatterplot(
    data = fov_meta, x = 'X_position', y = 'Y_position',
    hue = fov_cpm.to_df()[gene].values, linewidth = 0, s=5, palette = cmap,
    legend=None,
)
#%%
# %%
gene = 'CD21'
cutoff = 9.5
plt.hist(np.array(fov_cpm[:,gene].X.flatten()), bins=20)
fov_cpm.obs[gene + '_gate'] = (fov_cpm.to_df()[gene]>cutoff).values + 0
sc.pl.umap(
    fov_cpm, color = ['leiden',gene, gene + '_gate'], ncols=2, cmap='viridis', 
    vmin=[None, 'p5', None],
    vmax=[None, 'p95', None])

# %%
gene = 'CD20'
cutoff = 6
plt.hist(np.array(fov_cpm[:,gene].X.flatten()), bins=20)
fov_cpm.obs[gene + '_gate'] = (fov_cpm.to_df()[gene]>cutoff).values + 0
sc.pl.umap(
    fov_cpm, color = ['leiden',gene, gene + '_gate'], ncols=2, cmap='viridis', 
    vmin=[None, 'p5', None],
    vmax=[None, 'p95', None])
# %%
gene = 'CD14'
cutoff = 9.5
plt.hist(np.array(fov_cpm[:,gene].X.flatten()), bins=20)
fov_cpm.obs[gene + '_gate'] = (fov_cpm.to_df()[gene]>cutoff).values + 0
sc.pl.umap(
    fov_cpm, color = ['leiden',gene, gene + '_gate'], ncols=2, cmap='viridis', 
    vmin=[None, 'p5', None],
    vmax=[None, 'p95', None])
#%%
gene = 'CD45'
cutoff = 6.5
plt.hist(np.array(fov_cpm[:,gene].X.flatten()), bins=20)
fov_cpm.obs[gene + '_gate'] = (fov_cpm.to_df()[gene]>cutoff).values + 0
sc.pl.umap(fov_cpm, color = ['leiden',gene, gene + '_gate'], ncols=2, cmap='viridis')
# %%
gene = 'CD4'
cutoff = 9
plt.hist(np.array(fov_cpm[:,gene].X.flatten()), bins=20)
fov_cpm.obs[gene + '_gate'] = (fov_cpm.to_df()[gene]>cutoff).values + 0
sc.pl.umap(fov_cpm, color = ['leiden',gene, gene + '_gate'], ncols=2, cmap='viridis')
# %%
gene = 'CD8A'
cutoff = 8
plt.hist(np.array(fov_cpm[:,gene].X.flatten()), bins=20)
fov_cpm.obs[gene + '_gate'] = (fov_cpm.to_df()[gene]>cutoff).values + 0
sc.pl.umap(fov_cpm, color = ['leiden',gene, gene + '_gate'], ncols=2, cmap='viridis')
# %%
gene = 'CD3D'
cutoff = 7
plt.hist(np.array(fov_cpm[:,gene].X.flatten()), bins=20)
fov_cpm.obs[gene + '_gate'] = (fov_cpm.to_df()[gene]>cutoff).values + 0
sc.pl.umap(fov_cpm, color = ['leiden',gene, gene + '_gate'], ncols=2, cmap='viridis')
# %%
gene = 'CD68'
cutoff = 6.5
plt.hist(np.array(fov_cpm[:,gene].X.flatten()), bins=20)
fov_cpm.obs[gene + '_gate'] = (fov_cpm.to_df()[gene]>cutoff).values + 0
sc.pl.umap(fov_cpm, color = ['leiden',gene, gene + '_gate'], ncols=2, cmap='viridis')
#%%
gene = 'CD163'
cutoff = 9
plt.hist(np.array(fov_cpm[:,gene].X.flatten()), bins=20)
fov_cpm.obs[gene + '_gate'] = (fov_cpm.to_df()[gene]>cutoff).values + 0
sc.pl.umap(fov_cpm, color = ['leiden',gene, gene + '_gate'], ncols=2, cmap='viridis')
#%%
gene = 'CD11B'
cutoff = 9
plt.hist(np.array(fov_cpm[:,gene].X.flatten()), bins=20)
fov_cpm.obs[gene + '_gate'] = (fov_cpm.to_df()[gene]>cutoff).values + 0
sc.pl.umap(fov_cpm, color = ['leiden',gene, gene + '_gate'], ncols=2, cmap='viridis')
