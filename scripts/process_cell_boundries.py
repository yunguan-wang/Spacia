#%%
#%%
import pandas as pd
import numpy as np
import os
#%%
# Cell makss
input_path = '/project/shared/xiao_wang/projects/cell2cell_inter/data/merscope_data/HumanLungCancerPatient1'
cell_coords = pd.Series()
import h5py
boundries_fn = os.listdir(input_path + '/cell_boundaries')
for bfn in boundries_fn:
    print(bfn)
    bfn = os.path.join(input_path, 'cell_boundaries', bfn)
    f = h5py.File(bfn,'r')
    for cell in list(f['featuredata']):
        coords = np.array((f['featuredata'][cell]['zIndex_3']['p_0']['coordinates'][0]))
        if coords.shape[0] >= 5:  
            cell_coords[cell] = coords
    f.close()
cell_coords = cell_coords.to_frame(name='coord')
cell_coords['X'] = cell_coords.coord.apply(lambda x: '_'.join(x[:,0].round(2).astype(str)))
cell_coords['Y'] = cell_coords.coord.apply(lambda x: '_'.join(x[:,1].round(2).astype(str)))
cell_coords.iloc[:,1:].to_csv(input_path + '/cell_coords.csv')