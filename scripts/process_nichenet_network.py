#%%
from sys import meta_path
from matplotlib.style import available
import pandas as pd
import os
import matplotlib.pyplot as plt
# %%
def get_job_interactions(
    nichenet, avg_cpm, m_genes, s_genes, cell_pairs):
    '''
    get a table of interactions with data necessary for constructing spacia jobs.
    '''
    job_net = nichenet.copy()
    m_genes = [x for x in m_genes if x in avg_cpm.columns]
    s_genes = [x for x in s_genes if x in avg_cpm.columns]
    print('M:', len(m_genes), 'S:', len(s_genes))
    job_net = job_net[job_net.level_0.isin(avg_cpm.columns)]
    job_net = job_net[job_net.level_1.isin(avg_cpm.columns)]
    job_net.columns = ['RG','SG','Nichenet_weight']
    job_net['SG_type'] = ['M' if x in m_genes else 'S' for x in job_net.SG]
    job_interactions = pd.DataFrame()
    for cp in cell_pairs:
        ct1, ct2 = cp
        job_rgs = []
        job_sgs = []
        n_sg_used = 0
        for sg_list in [m_genes,s_genes]:
            for sg in sg_list:
                if n_sg_used == 20:
                    break
                if avg_cpm.loc[ct1, sg]<avg_cpm[sg].median():
                    continue
                net = job_net[job_net.SG==sg].copy()
                net['Sender'] = ct1
                net['Receiver'] = ct2
                net['cp'] = ct1 + '_' + ct2
                rg = net.RG.values
                rg = rg[avg_cpm.loc[ct2, rg]> avg_cpm[rg].median()][:5]
                job_interactions = job_interactions.append(
                    net[net.RG.isin(rg)])
                job_sgs.append(sg)
                job_rgs += rg.tolist()
                n_sg_used+=1
    job_interactions = job_interactions.reset_index().iloc[:,1:]
    return job_interactions

def job_writer(
    job_interactions, 
    output_fd,
    cpm_path,
    meta_path,
    spacia_path = '/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/code/scripts/spacia.py'
    ):
    '''
    Contruct spacia jobs for each cell type pair
    '''
    output_fd = os.path.join(output_fd, 'spacia')
    if not os.path.exists(output_fd):
        os.makedirs(output_fd)
    for _, g_df in job_interactions.groupby('cp'):
        job_name = (g_df.Sender + '_vs_' + g_df.Receiver).iloc[0].replace('/','OR')
        job_path = os.path.join(output_fd, job_name)
        sc, rc = g_df.iloc[0,4:6]
        sf = ','.join(g_df.SG.unique())
        rf = ','.join(g_df.RG.unique())
        with open(os.path.join(output_fd, '{}.sh'.format(job_name)), 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#SBATCH --job-name=spacia \n')
            f.write('#SBATCH --partition=256GB\n')
            f.write('#SBATCH --nodes=1\n')
            f.write('#SBATCH --cpus-per-task=32\n')
            f.write('#SBATCH --time=10-00:00:00\n')
            f.write('#SBATCH --output=./{}.log\n'.format(job_name))
            f.write('#SBATCH --error=./{}.err\n'.format(job_name))
            
            f.write('source /home2/s190548/.bashrc\n')
            f.write('conda activate p3s\n')
            f.write('module load R/4.0.2-gccmkl\n')
            k=0
            f.write(
                'python {} {} {} -rc {} -sc {} -sf {} -rf {} -nc 20 -o {}\n'.format(
                    spacia_path, cpm_path, meta_path, rc, sc, sf, rf, job_path
                    )
                )
# %%
input_path = '/project/shared/xiao_wang/projects/cell2cell_inter/data/NicheNet'
output_path = '/project/shared/xiao_wang/projects/cell2cell_inter/data/spacia_merscope/nichenet_test'
os.chdir(input_path)
# %%
# Process nichenet
nichenet = pd.read_csv(
    'NicheNet_ligand_target.txt', sep='\t', index_col=0
)
ignored_genes = [x for x in nichenet.index if x[:4] == 'LINC']
ignored_genes += [x for x in nichenet.index if x[:3] == 'MIR']
ignored_genes += [x for x in nichenet.index if x[-4:] == '-AS1']
ignored_genes += [x for x in nichenet.index if x[:4] == 'LOC1']
nichenet = nichenet[~nichenet.index.isin(ignored_genes)]
top_net = nichenet.stack().sort_values(ascending=False)[:50000]
top_net = top_net.reset_index()
#%%
cell_pairs = [
    ('Granulocytes','Tumor_epithelial_cells'),
    ('Macrophages/Dendritic_cells','Tumor_epithelial_cells'),
    ('Fibroblasts','Tumor_epithelial_cells'),
    ('CD8_T_cells','Tumor_epithelial_cells'),
    ('Endothelial_cells','Tumor_epithelial_cells'),
    ('NK_cells','Tumor_epithelial_cells'),
    ('Th_cells','Tumor_epithelial_cells'),
    ('B_cells','Tumor_epithelial_cells'),
    ('Granulocytes','Normal_epithelial_cells'),
    ('Macrophages/Dendritic_cells','Normal_epithelial_cells'),
    ('Fibroblasts','Normal_epithelial_cells'),
    ('CD8_T_cells','Normal_epithelial_cells'),
    ('Endothelial_cells','Normal_epithelial_cells'),
    ('NK_cells','Normal_epithelial_cells'),
    ('Th_cells','Normal_epithelial_cells'),
    ('B_cells','Normal_epithelial_cells'),
    ('Th_cells','CD8_T_cells'),
    ('Th_cells', 'B_cells'),
    ('Treg_cells', 'CD8_T_cells'),
    ('Treg_cells', 'Th_cells'),
]
top_ligands = nichenet.mean().sort_values(ascending=False).index[:200]
original_order = top_ligands
top_ligands = pd.DataFrame(index = top_ligands)
uniprot = pd.read_csv('/project/shared/xiao_wang/projects/cell2cell_inter/data/uniprot_protein_metadata.tsv', sep='\t', index_col=2)
top_ligands = top_ligands.merge(
    uniprot[['Subcellular location [CC]']],
    left_index=True, right_index=True, sort=False, how='left')
top_ligands['Secreted'] = top_ligands['Subcellular location [CC]'].str.contains('Secreted')
top_ligands['Secreted'] = top_ligands['Secreted'].fillna(True)
top_ligands = top_ligands.loc[original_order]
secreted = top_ligands[top_ligands.Secreted].index
membrane = top_ligands[~top_ligands.Secreted].index

merscope_path = '/project/shared/xiao_wang/projects/cell2cell_inter/data/merscope_data/'
merscope_dsets = [merscope_path + x for x in os.listdir(merscope_path)]
cell_types = list(set([x for sublist in cell_pairs for x in sublist]))
for fd in merscope_dsets:
    tmp_path = os.path.join(merscope_path, fd, 'spacia_spot_meta.txt')
    if not os.path.exists(tmp_path):
        print('{} do not exist!'.format(tmp_path))
        continue
    meta_path = tmp_path
    cpm_path = os.path.join(merscope_path, fd, 'sprod/denoised_stiched.txt')
    cpm = pd.read_csv(cpm_path, index_col=0, sep='\t')
    spot_meta = pd.read_csv(meta_path, index_col=0, sep='\t')
    avg_cpm = cpm.groupby(spot_meta.cell_type).mean()
    un_matched = [x for x in cell_types if x not in avg_cpm.index]
    if len(un_matched)>0:
        print(un_matched)
        break
    job_interactions = get_job_interactions(
        top_net, avg_cpm, membrane, secreted, cell_pairs)
    job_interactions.to_csv(
        os.path.join(fd,'spacia','job_nichenet_interactions.csv'))
    job_writer(job_interactions, fd, cpm_path, meta_path)

# %%


#%%
# def job_writer(job_f, pairs, cp):
#     '''
#     each pair is in the order of sender, receiver
#     '''
#     sc, rc = cp
#     spacia_path = '/endosome/work/InternalMedicine/s190548/software/cell2cell_inter/code/scripts/spacia.py'
#     with open(os.path.join(output_path, 'job_spacia.sh'), 'w') as f:
#         f.write('#!/bin/bash\n')
#         f.write('#SBATCH --job-name=sprod \n')
#         f.write('#SBATCH --partition=256GB\n')
#         f.write('#SBATCH --nodes=1\n')
#         f.write('#SBATCH --cpus-per-task=32\n')
#         f.write('#SBATCH --time=10-00:00:00\n')
#         f.write('#SBATCH --output=./spacia_nichenet.log\n')
#         f.write('#SBATCH --error=./spacia_nichenet.err\n')
        
#         f.write('source /home2/s190548/.bashrc\n')
#         f.write('conda activate p3s\n')
#         f.write('module load R/4.0.2-gccmkl\n')
#         k=0
#         for p in pairs:
#             job_path = os.path.join(output_path, '_'.join(p))
#             if k == 10:
#                 sep = '\n'
#                 k=0
#             else:
#                 sep = '&\n'
#             f.write(
#                 'python {} {} {} -rc {} -sc {} -sf {} -rf {} -nc 20 -o {} {}'.format(
#                     spacia_path, cpm_path, meta_path,rc, sc, p[0], p[1], job_path, sep
#                     )
#                 )
#             k+=1
# %%
