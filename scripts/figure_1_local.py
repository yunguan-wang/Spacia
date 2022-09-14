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

def marker_overlay(
    gene1, gene2, spot_meta, cpm, img, cells_1, cells_2, color1, color2, interactions, 
    output=None, no_raw_img = False):
    spot_meta = spot_meta.copy()
    x_min, x_max = spot_meta.X.min()//100 * 100, (1+spot_meta.X.max()//100) * 100
    y_min, y_max = spot_meta.Y.min()//100 * 100, (1+spot_meta.Y.max()//100) * 100
    gene_exp_mask = np.zeros(img.shape, dtype = 'uint8')
    gene_exp_mask[:,:,:] = 255
    r = spot_meta['Spot_radius'][0]
    bbox_mask = morphology.disk(r)

    for i in range(2):
        if i == 0:
            cells, gene, color = cells_1, gene1, color1
        else:
            cells, gene, color = cells_2, gene2, color2
        gene_v = cpm.loc[cells, gene]
        gene_v = (255*gene_v / gene_v.max()).astype('uint8')
        gene_v[gene_v == 0] = 20
        cmap = sns.light_palette(color, as_cmap=True)
        norm = mpl.colors.Normalize(vmin=0, vmax=255)
        mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
        for id, row in spot_meta.iterrows():
            x, y = row[['X','Y']].astype(int)
            if id in cells:
                gene_exp_level = gene_v[id]
                rgb_v = (255*np.array(mapper.to_rgba(gene_exp_level))[:3]).astype('uint8')
                gene_exp_mask[x-r:x+r+1,y-r:y+r+1, :][bbox_mask==1] = rgb_v
    gene_exp_mask = gene_exp_mask.astype('uint8')
    # tmp_array = np.zeros(img.shape, 'uint8')
    # tmp_array[:,:,1] = gene_exp_mask
    # gene_exp_mask = tmp_array
    _, ax = plt.subplots(figsize=(36,36))
    if no_raw_img:
        merged_img = gene_exp_mask
        io.imshow(merged_img[:,:,1], cmap="gray", ax = ax)
    else:
        merged_img = (
            0.3*img[x_min:x_max, y_min:y_max] + 0.7*gene_exp_mask[x_min:x_max, y_min:y_max]
            ).astype('uint8')
        io.imshow(merged_img, ax = ax)
    spot_meta['X'] = spot_meta['X'] - x_min
    spot_meta['Y'] = spot_meta['Y'] - y_min
    for _, row in interactions.iterrows():
        try:
            x, y = spot_meta.loc[row['sender'], ['X','Y']]
            x1, y1 = spot_meta.loc[row['receiver'], ['X','Y']]
            ax.arrow(y,x, y1-y, x1-x, width=1, color='k')
        except:
            continue
    plt.axis('off')
    if output is not None:
        plt.savefig(output)
        plt.close()
    return gene_exp_mask

def cal_qc(adata):
    adata.var_names_make_unique()
    # add the total counts per cell as observations-annotation to adata
    adata.obs["n_counts"] = adata.X.sum(axis=1)
    adata.obs["n_genes"] = np.sum(adata.X>0,axis=1)
    return adata

def normalize_adata(adata):
    # sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e3)
    sc.pp.log1p(adata)
    adata.raw = adata
    return adata

def gm_gating(df, nc = 3):
    gm = GaussianMixture(n_components=nc, init_params='kmeans', random_state=0)
    labels = gm.fit_predict(df)
    target_cluster = (labels == df.groupby(labels).mean().idxmax()[0])
    return df.index[target_cluster]

def percentile_gating(df, gene, q_cutoff=0.5):
    cutoff = df[gene].quantile(q_cutoff)
    pop = df[gene].groupby(cts_sub.obs.leiden.astype(str)).median() > cutoff
    pop_cells = pop.index[pop.values]
    return pop, pop_cells

def plot_product_dist(
    df, gene1, gene2, meta_data, x_col = 'center_x', y_col = 'center_y'):
    keep = (df[[gene1, gene2]]>0).any(axis=1)
    df = df.loc[keep, [gene1, gene2]]
    df_meta = meta_data.loc[df.index, [x_col, y_col]]
    gene_product = np.matmul(df[[gene1]].values, df[[gene2]].T.values)
    np.fill_diagonal(gene_product, np.nan)
    dists = cdist(df_meta, df_meta)
    np.fill_diagonal(dists, np.nan)
    product_col = gene1 + '_' + gene2 + '_product'
    plot_data = pd.DataFrame(gene_product.flatten(), columns = [product_col])
    plot_data['distance'] = dists.flatten()
    plot_data = plot_data.dropna()
    # distance_bins = pd.cut(
    #     plot_data.distance, bins=10,
    #     labels = ['Distance_bin_' + str(i) for i in range(1,11)])
    # plot_data['Distance_bins'] = distance_bins.copy()
    plot_data['Product_bins'] = pd.cut(
        plot_data[product_col], bins=10,
        labels = ['Product_bin_' + str(i) for i in range(1,11)])
    plot_data[product_col] = plot_data[product_col].astype(float)
    # plot_data = plot_data[plot_data[product_col]>0]
    # sns.boxplot(data = plot_data, x = 'Distance_bins', y = product_col,)
    # plt.xticks(rotation=90)
    sns.boxplot(data = plot_data, x = 'Product_bins', y = 'distance',)
    plt.xticks(rotation=90)

def plot_product_dist_by_type(
    df1, df2, meta_data, x_col = 'center_x', y_col = 'center_y'):
    gene_product = np.matmul(df1.values.reshape(-1,1), df2.values.reshape(1,-1))
    coords = meta_data[[x_col, y_col]]
    dists = cdist(coords.loc[df1.index], coords.loc[df2.index])
    product_col = df1.name + '_' + df2.name + '_product'
    plot_data = pd.DataFrame(gene_product.flatten(), columns = [product_col])
    plot_data['distance'] = dists.flatten()
    distance_bins = pd.cut(
        plot_data.distance, bins=[-0.1,100, 500, 1000, 5000],
        labels = ['Distance_bin_' + str(i) for i in range(1,5)])
    plot_data['Distance_bins'] = distance_bins.copy()
    # plot_data = plot_data[plot_data['distance']<2000]
    # plot_data['Product_bins'] = pd.cut(
    #     plot_data[product_col], bins=[-0.1, 1, 5, 9, 15],
    #     labels = ['Product_bin_' + str(i) for i in range(1,5)])
    # plot_data['Product_bins'] = pd.cut(
    #     plot_data[product_col], bins=[-0.1,2,6,10,14],
    #     labels = ['Product_bin_' + str(i) for i in range(1,5)])
    plot_data[product_col] = plot_data[product_col].astype(float)
    # plot_data = plot_data[plot_data[product_col]>0]
    sns.boxplot(data = plot_data, x = 'Distance_bins', y = product_col,)
    plt.xticks(rotation=90)
    # sns.boxplot(data = plot_data, x = 'Product_bins', y = 'distance',)
    # sns.violinplot(
    #     data = plot_data, x = 'Product_bins', y = 'distance',cut=0, 
    #     inner='quartile')
    # plt.xticks(rotation=90)

def plot_pooled_product_dist_by_type(
    sender,
    receiver,
    meta_data,
    x_col = 'center_x',
    y_col = 'center_y',
    distcutoffs = [0, 10, 20, 40, 80, 200],
    receiver_gene_qcuts = [0, 0.5, 1],
    ):
    gene1 = sender.name
    gene2 = receiver.name
    bin2 = receiver.quantile(receiver_gene_qcuts)
    coords = meta_data[[x_col, y_col]]
    dists = cdist(coords.loc[sender.index], coords.loc[receiver.index])
    dists = pd.DataFrame(dists, index = sender.index, columns = receiver.index)
    plot_data = dists.stack().reset_index()
    plot_data.columns = ['sender', 'receiver', 'distance']
    plot_data = plot_data.dropna()
    plot_data[gene1] = sender.loc[plot_data.sender].values
    plot_data[gene2] = receiver.loc[plot_data.receiver].values
    plot_data['Distance_bin'] = pd.cut(
        plot_data.distance, bins=distcutoffs, include_lowest = True)
    plot_data[gene2 + '_bin'] = pd.cut(
        plot_data[gene2], bins = bin2, include_lowest = True, duplicates='drop')
    plot_data = plot_data.dropna()
    gene1sum = plot_data.groupby(
        ['Distance_bin', 'receiver']
        )[gene1].mean().reset_index()
    gene1sum[gene1] = gene1sum[gene1].fillna(0)
    gene1sum.columns = ['Distance_bin', 'receiver'] + [gene1 + '_sum']
    plot_data = plot_data.merge(gene1sum, on=['Distance_bin', 'receiver'])
    plot_data['Nearby_'+ gene1 + '_expression'] = pd.qcut(
        plot_data[gene1 + '_sum'],q = 2, duplicates='drop')
    sns.boxplot(
        data = plot_data,
        x = 'Distance_bin',
        y = gene2,
        hue = 'Nearby_'+ gene1 + '_expression')
    plt.legend(
        bbox_to_anchor = (1, 0.5), loc = 'center left',
        title= 'Nearby_'+ gene1 + '_expression')
    # sns.boxplot(
    #     data = plot_data,
    #     x = 'Distance_bin',
    #     y = gene1 + '_sum',
    #     hue = gene2 + '_bin')
    # plt.legend(
    #     bbox_to_anchor = (1, 0.5), loc = 'center left',
    #     title= gene2 + '_bin')
    plt.xticks(rotation=90)
    return plot_data

def corr_product_dist(
    df, gene1, gene2, meta_data, x_col = 'center_x', y_col = 'center_y'):
    keep = (df[[gene1, gene2]]>0).any(axis=1)
    df = df.loc[keep, [gene1, gene2]]
    df_meta = meta_data.loc[df.index, [x_col, y_col]]
    gene_product = np.matmul(df[[gene1]].values, df[[gene2]].T.values)
    np.fill_diagonal(gene_product, np.nan)
    dists = cdist(df_meta, df_meta)
    np.fill_diagonal(dists, np.nan)
    dist_vector = dists.flatten()[~np.isnan(dists.flatten())]
    product_vector = gene_product.flatten()[~np.isnan(gene_product.flatten())]
    return gene1, gene2, np.corrcoef(dist_vector, product_vector)[0,1]
#%%    
input_path = '/project/shared/xiao_wang/projects/cell2cell_inter/data/merscope_data/HumanLungCancerPatient1'
output_path = '/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/figure1'
sc.settings.figdir = output_path
sns.set_theme(
    context='paper', style='white', palette='gist_rainbow', 
    font= 'Arial',font_scale=3)
sns.set_style({'ytick.left': True, 'xtick.bottom': True})
#%%
cpm = pd.read_csv(
    input_path + '/sprod/denoised_stiched.txt', index_col=0, sep='\t')
# cpm = pd.read_csv(
#     input_path + '/sprod/Counts.txt', index_col=0, sep='\t')
meta_data = pd.read_csv(input_path + '/sprod/Spot_metadata.csv', index_col=0)
cell_types = pd.read_csv(
    input_path + '/typed_cells_noah.txt', index_col=0, sep= '\t', header=None
    ).index
cell_types = [
    'Myeloid cells' if x in [
        'Granulocytes', 'Macrophages/Dendritic cells'
        ] else x for x in cell_types
        ]
umap_x = pd.read_csv(input_path + '/umap_embedding.csv')
umap_x.index = ['cell_' + str(i+1) for i in umap_x.index]
meta_data.loc[umap_x.index, 'celltype'] = cell_types

#%%
gene1 = 'TGFB'
max_d = 30
fib = meta_data[meta_data.celltype == 'Fibroblasts']
tumor = meta_data[meta_data.celltype == 'Tumor epithelial cells']
fib_to_tumor_dist = cdist(fib[['X','Y']], tumor[['X','Y']])
caf = fib.index[(fib_to_tumor_dist<=100).any(axis=1)]
receiver_cells = tumor.index
sender = cpm.loc[caf,['TGFB1', 'TGFB2','TGFB3']].mean(axis=1)
sender.name = 'TGFB'
coords = meta_data[['X', 'Y']]

dists = cdist(coords.loc[sender.index], coords.loc[receiver_cells])
crit1 = (dists <= max_d).any(axis=0)
crit2 = (dists <= max_d).any(axis=1)
dists = dists[crit2,:][:,crit1]
dists = pd.DataFrame(
    dists, index = sender.index[crit2], columns = receiver_cells[crit1])
dists = dists.stack().reset_index()
dists.columns = ['sender', 'receiver', 'distance']
dists = dists[dists.distance<=max_d]
dists['Distance_bin'] = pd.cut(
    dists.distance, bins=[0,7,max_d], include_lowest = True)
dists = dists.dropna()

emt_genes = [
    'CDH2','CDH11','FN1','VIM','TWIST1','SNAI1','ZEB1','ZEB2','DCN',
    'SNAI2','DSP', 'MMP2', 'MMP9']
emt_genes = [x for x in emt_genes if x in cpm.columns]

for gene2 in emt_genes:
    plot_data = dists.copy()
    plot_data[gene1] = sender.loc[plot_data.sender.values].values
    plot_data[gene2] = cpm.loc[plot_data.receiver.values, gene2].values

    gene1sum = plot_data.groupby(
        ['Distance_bin', 'receiver']
        )[gene1].mean().reset_index()
    gene1sum[gene1] = gene1sum[gene1].fillna(0)
    gene1sum.columns = ['Distance_bin', 'receiver'] + [gene1 + '_sum']
    plot_data = plot_data.merge(gene1sum, on=['Distance_bin', 'receiver'])
    plot_data = plot_data.drop_duplicates(subset = ['Distance_bin','receiver'])
    sender_gene_col = 'Fibroblast_' + gene1 + '_expression'
    for _ ,df in plot_data.groupby('Distance_bin'):
        cuts = pd.qcut(
            df[gene1 + '_sum'],q = 2, duplicates='drop',
            labels = ['Low','High'])
        plot_data.loc[
            cuts.index, sender_gene_col
            ] = cuts.values.astype(str)

    cmap = ['r','cyan']
    _ = plt.figure(figsize=(10,6))
    ax = sns.boxplot(
        data = plot_data,
        x = 'Distance_bin',
        y = gene2,
        hue = sender_gene_col,
        hue_order=['Low','High'],
        palette = cmap,
        linewidth=5)
    plt.legend(
        bbox_to_anchor = (1, 0.5), loc = 'center left',
        title= sender_gene_col.replace('_','\n'))
    plt.ylabel('Tumor ' + gene2 + ' expression', fontdict={'weight':'bold'})
    plt.xlabel('')
    plt.xticks(
        [0,1], labels = ['Adjacent', 'Distant'], fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path + '/CAF_tumor_TGFB_' + gene2 + '_raw.pdf')
    plt.close()
# %%
pseudo_img = meta_data.sample(50000, replace=False, random_state=0)
pseudo_img = pd.concat([pseudo_img,cpm], axis=1, join = 'inner')
_ = plt.figure(figsize=(16,10))
colors = sns.color_palette(palette='hsv', n_colors=10, as_cmap=False).as_hex()
colors.append('#808080')
color_order = [
    'Tumor epithelial cells',
    'Normal epithelial cells',
    'CD8 T cells',
    'Th cells',
    'Treg cells',
    'NK cells',
    'B cells',
    'Endothelial cells',
    'Fibroblasts',
    'Myeloid cells',
    'unknown',
    ]
sns.scatterplot(
    data = pseudo_img, x = 'X', y = 'Y', hue = 'celltype', linewidth=0, 
    s = 10, palette = colors, hue_order = color_order)
plt.legend(bbox_to_anchor = (1, 0.5), loc='center left', markerscale=3)
plt.tight_layout()
plt.savefig(output_path + '/celltypes_pseudoimage.pdf')

# %%
# ======================================================================
# ======================================================================
# ======================================================================
# Visium
input_path = '/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/'
output_path = '/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/figure1'
# input_path = 'P:/software/cell2cell_inter/data/'
# output_path = 'P:/software/cell2cell_inter/data/figure1'
sc.settings.figdir = output_path
sns.set_theme(
    context='paper', style='white', palette='gist_rainbow', 
    font= 'Arial',font_scale=3)
sns.set_style({'ytick.left': True, 'xtick.bottom': True})
# %%
cpm = pd.read_csv(
    input_path + '/visium_denoised.txt', index_col=0, sep='\t')
cpm = np.log1p(cpm.apply(lambda x: 1e4 * x/x.sum()))
# cpm = pd.read_csv(
#     input_path + '/sprod/Counts.txt', index_col=0, sep='\t')
meta_data = pd.read_csv(input_path + '/visium_celltyping_metadata.csv', index_col=0)
#%%
gene1 = 'PDCD1'
max_d = 10000
stromal = meta_data[meta_data['Type'] == 'I/S']
tumor = meta_data[meta_data['Type'] != 'I/S']
t2s_dist = cdist(stromal[['X','Y']], tumor[['X','Y']])
stromal = stromal.index[(t2s_dist<=max_d).any(axis=1)]
receiver_cells = tumor.index
sender_cells = stromal
sender_exp = cpm.loc[sender_cells, gene1]
coords = meta_data[['X', 'Y']]

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
    dists.distance, bins=[0,280,600, max_d], include_lowest = True)
dists = dists.dropna()

#%%
gene2 = 'CD274'
plot_data = dists.copy()
plot_data[gene1] = cpm.loc[plot_data.sender.values, gene1].values
plot_data[gene2] = cpm.loc[plot_data.receiver.values, gene2].values

gene1sum = plot_data.groupby(
    ['Distance_bin', 'receiver']
    )[gene1].mean().reset_index()
gene1sum[gene1] = gene1sum[gene1].fillna(0)
gene1sum.columns = ['Distance_bin', 'receiver'] + [gene1 + '_sum']
plot_data = plot_data.merge(gene1sum, on=['Distance_bin', 'receiver'])
plot_data = plot_data.drop_duplicates(subset = ['Distance_bin','receiver'])
sender_gene_col = gene1 + '_sum'
for _ ,df in plot_data.groupby('Distance_bin'):
    cuts = pd.qcut(
        df[gene1 + '_sum'],q = 2, duplicates='drop',
        labels = ['Low','High'])
    plot_data.loc[
        cuts.index, sender_gene_col
        ] = cuts.values.astype(str)

cmap = ['r','cyan']
_ = plt.figure(figsize=(10,6))
ax = sns.boxplot(
    data = plot_data,
    x = 'Distance_bin',
    y = gene2,
    hue = sender_gene_col,
    hue_order=['Low','High'],
    palette = cmap,
    linewidth=5)
plt.legend(
    bbox_to_anchor = (1, 0.5), loc = 'center left',
    title= 'Nearby Stromal\nPDCD1 expression')
plt.ylabel('Tumor ' + gene2 + ' expression', fontdict={'weight':'bold'})
plt.xlabel('')
plt.xticks(
    [0,1], labels = ['Adjacent', 'Distant'], fontweight='bold')
plt.tight_layout()
# plt.savefig(output_path + '/Stromal_tumor_' + gene1 + '_' + gene2 + '.pdf')
# plt.close()
# %%
for gene in ['CDH1', 'FN1', 'PTPRC','celltype']:
    _ = plt.figure(figsize=(14,10))
    sns.scatterplot(
        data = pseudo_img, x = 'X', y = 'Y', hue = gene, linewidth=0, 
        s = 10)
    plt.legend(bbox_to_anchor = (1, 0.5), loc='center left', markerscale=3)
    plt.title(gene)
    plt.tight_layout()
    plt.savefig(output_path + '/' + gene + '_pseudoimage.png')
    plt.close()
    
#%%
meta_data = pd.read_csv(input_path + '/visium_celltyping_metadata.csv', index_col=0)
img = io.imread(
    '/project/shared/xiao_wang/projects/MOCCA/data/Visium/breast_cancer/V1_Breast_Cancer_Block_A_Section_1_image.tif')

gene2 = 'CD274'
interactions = dists.copy()
interactions[gene1] = cpm.loc[interactions.sender.values, gene1].values
interactions[gene2] = cpm.loc[interactions.receiver.values, gene2].values

gene1sum = interactions.groupby(
    ['Distance_bin', 'receiver']
    )[gene1].mean().reset_index()
gene1sum[gene1] = gene1sum[gene1].fillna(0)
gene1sum.columns = ['Distance_bin', 'receiver'] + [gene1 + '_sum']
interactions = interactions.merge(gene1sum, on=['Distance_bin', 'receiver'])
cells1 = interactions[interactions.distance <= 600].sender.unique()
cells2 =  interactions[interactions.distance <= 600].receiver.unique()
interactions = interactions[interactions.distance<=280]
plot_meta = meta_data[
    (meta_data.X>=15000) &
    (meta_data.X<=21000) &
    (meta_data.Y>=10000) &
    (meta_data.Y<=19000) 
    ].copy()
# pdcd1_sum = gene1sum[gene1sum['Distance_bin'] == gene1sum['Distance_bin'].unique()[0]]
# cpm.loc[pdcd1_sum.receiver.values, 'PDCD1_sum'] = pdcd1_sum.PDCD1_sum.values
# # cells1 = plot_meta[plot_meta['Type'] == 'I/S'].index
# # cells2 =  plot_meta[plot_meta['Type'] != 'I/S'].index
# _ = marker_overlay(
#     'PDCD1', 'PDCD1_sum',
#     plot_meta, cpm, img=img, 
#     cells_1=cells1,
#     cells_2=cells2,
#     interactions=interactions,
#     output=output_path + '/visium_figure_PDCD1_PDCD1_sum_zoom.pdf',
#     )

_ = marker_overlay(
    'PDCD1', 'CD274',
    plot_meta, cpm, img=img, 
    cells_1=cells1,
    cells_2=cells2,
    color1='#79af97',
    color2='#df8f44',
    interactions=interactions,
    output=output_path + '/visium_figure_PDCD1_CD274_zoom.pdf',
    )

_ = marker_overlay(
    'PDCD1', 'CD274',
    meta_data, cpm, img=img, 
    cells_1=cells1,
    cells_2=cells2,
    color1='#79af97',
    color2='#df8f44',
    interactions=interactions,
    output=output_path + '/visium_figure_PDCD1_CD274.pdf',
    )

