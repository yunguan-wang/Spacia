import pandas as pd
import numpy as np
import os
#%%
# input_path = '/project/shared/xiao_wang/projects/cell2cell_inter/data/merscope_data/HumanLungCancerPatient1'
# output_path = '/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/data/figure1'

# # prepare spacia inputs
# meta_data = pd.read_csv(input_path + '/cell_metadata.csv', index_col=0)
# cell_types = pd.read_csv(
#     input_path + '/typed_cells_noah.txt', index_col=0, sep= '\t', header=None).index
# meta_data['cluster'] = cell_types.values
# spacia_meta_data = meta_data.iloc[:,[2,3,-1]]
# spacia_meta_data.columns = ['X','Y','cluster']
# spacia_meta_data.to_csv(input_path + '/spacia_spot_meta.txt', sep='\t')

# raw_cts = pd.read_csv(input_path + '/cell_by_gene.csv', index_col=0)
# blank_genes = [x for x in raw_cts.columns if x[:5] == 'Blank']
# raw_cts = raw_cts.drop(blank_genes, axis=1)
# raw_cts.index = ['cell_' + str(i+1) for i in raw_cts.index]
# umap_x = pd.read_csv(input_path + '/umap_embedding.csv')
# umap_x.index = ['cell_' + str(i+1) for i in umap_x.index]
# raw_cts = raw_cts.loc[umap_x.index]
# raw_cts.to_csv(input_path, '/spacia_counts.txt', sep='\t')

#%%
input_path = '/project/shared/xiao_wang/projects/cell2cell_inter/data/merscope_data/'
folders = os.listdir(input_path)
sprod_path = '/project/shared/xiao_wang/projects/MOCCA/code/Spatial_denoise/sprod.py'
for fd in folders:
    fd = os.path.join(input_path, fd)
    fns = os.listdir(fd)
    output_cts_fn = os.path.join(fd, 'sprod', 'Counts.txt')
    output_meta_fn = os.path.join(fd, 'sprod', 'Spot_metadata.csv')
    if 'sprod' in fns:
        continue
    elif (
        ('cell_by_gene.csv' not in fns) | 
        ('umap_embedding.csv' not in fns) |
        ('cell_metadata.csv' not in fns)
        ):
        continue
    output_path = os.path.join(fd, 'sprod')
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    if os.path.exists(output_cts_fn):
        continue 
    print(fd)

    umap_x = pd.read_csv(fd + '/umap_embedding.csv')
    umap_x.index = ['cell_' + str(i+1) for i in umap_x.index]

    meta_data = pd.read_csv(fd + '/cell_metadata.csv', index_col=0)
    meta_data.index = ['cell_' + str(i+1) for i in meta_data.index]
    meta_data = meta_data.iloc[:,2:4]
    meta_data.columns = ['X','Y']
    meta_data = meta_data.loc[umap_x.index]

    raw_cts = pd.read_csv(fd + '/cell_by_gene.csv', index_col=0)
    blank_genes = [x for x in raw_cts.columns if x[:5] == 'Blank']
    raw_cts = raw_cts.drop(blank_genes, axis=1)
    raw_cts.index = ['cell_' + str(i+1) for i in raw_cts.index]
    raw_cts = raw_cts.loc[umap_x.index]

    raw_cts = raw_cts.apply(lambda x: 10000*x/x.sum(), axis=1)
    cpm = np.log2(1+raw_cts)
    
    # write job file
    with open(os.path.join(fd, 'sprod','run_sprod.sh'), 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('#SBATCH --job-name=sprod \n')
        f.write('#SBATCH --partition=256GB\n')
        f.write('#SBATCH --nodes=1\n')
        f.write('#SBATCH --cpus-per-task=32\n')
        f.write('#SBATCH --time=10-00:00:00\n')
        f.write('#SBATCH --output=./sprod_%j.log\n')
        f.write('#SBATCH --error=./sprod_%j.err\n')
        f.write('source /home2/s190548/.bashrc\n')
        f.write('conda activate mocca\n')
        f.write('module load R/4.0.2-gccmkl\n')
        f.write(
            'python {} {} {} -sn 50 -d -dg -pb 15 --input_type batch\n'.format(
            sprod_path,
            output_path,
            output_path,
            )
        )

    cpm.round(2).to_csv(output_cts_fn, sep='\t')
    meta_data.to_csv(output_meta_fn)
# %%
