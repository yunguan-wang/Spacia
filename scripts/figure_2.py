#%%
import pandas as pd
import scanpy as sc
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
# %%
input_path = '/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/figure1/merscope_lungcancer/spacia/'
output_path = '/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/figure1/merscope_lungcancer/spacia/'
sc.settings.figdir = output_path
sns.set_theme(
    context='paper', style='white', palette='gist_rainbow', 
    font= 'Arial',font_scale=3)
sns.set_style({'ytick.left': True, 'xtick.bottom': True})
# %%
interactions = pd.read_csv(input_path + 'Interactions.csv', index_col=0)
b_fdr = pd.read_csv(input_path + 'B_and_FDR.csv', index_col=0)
betas = pd.read_csv(input_path + 'Pathway_betas.csv', index_col=0)
receiver_modules = pd.read_json(
    os.path.join(input_path, 'model_input/receiver_pathways.json'), 
    orient='index')
sender_modules = pd.read_json(
    os.path.join(input_path, 'model_input/sender_pathways.json'), 
    orient='index')
# %%
sender_modules[sender_modules.apply(lambda x: 'TGFB2' in x.values, axis=1)]
# %%
df_int = pd.pivot(betas, columns = 'Sender_pathway', values = 'Beta')
df_int.index = [x.replace('module', 'r') for x in df_int.index]
df_int.columns = [x.replace('module', 's') for x in df_int.columns]
sns.clustermap(
    df_int, cmap='coolwarm', vmax=0.5, vmin=-0.5, figsize=(16,16), method='complete',
    xticklabels=1,
    yticklabels=1,
    annot_kws={"size": 2}
    )
# %%