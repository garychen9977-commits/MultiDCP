import pandas as pd
import h5py
import numpy as np
import pickle

data = pd.read_csv('Bayesian_GSE92742_with_confidence.csv')
# ## random_sel_dis for each sig_id
# data['rand_distil_id'] = data['distil_id'].str.split('|').str[0].str.replace(':', '_')

with h5py.File('Bayesian_GSE92742_Level5_COMPZ_n178984x978.h5', 'r') as level5_file:
    level5_df = pd.DataFrame(np.array(level5_file['data']), index = np.array(level5_file['colid']), columns = np.array(level5_file['rowid']))

## select qualified data
data = data[data['sig_id'].str.contains('24H')]
data = data[data['pert_type'] == 'trt_cp']

all_drugs = set(pd.read_csv('../DeepCE/DeepCE/data/all_drugs_l1000.csv').broad_cpd_id)
data = data[data['pert_id'].isin(all_drugs)]

all_cells = set(pd.read_csv('../DeepCE/DeepCE/data/adjusted_ccle_tcga_ad_tpm_log2.csv', index_col = 0).index)
print('{0!r} cells are not found in CCLE'.format(len(set(data.cell_id)-all_cells)))
data = data[data['cell_id'].isin(all_cells)]

all_doses = {"0.04 um", "0.12 um", "0.37 um", "1.11 um", "3.33 um", "10.0 um"}
data.pert_idose = data.pert_idose.str.replace('10 um', '10.0 um')
new_data = data[data['pert_idose'].isin(all_doses)]
print('{0!r} data was dropped because of the doses'.format(len(data)-len(new_data)))

merged_data = new_data.merge(level5_df, left_on = 'sig_id', right_index = True)
cols = level5_df.columns
merged_data = merged_data[['sig_id','pert_id', 'pert_type', 'cell_id', 'pert_idose'] + list(cols)]

bench = pickle.load(open('./benchmark_sigid.p', 'rb'))
test_data = merged_data[~merged_data.sig_id.isin(bench)]

