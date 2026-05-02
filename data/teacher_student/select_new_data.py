import pandas as pd
from scipy.stats import pearsonr
import os
import numpy as np
import sys
import pdb
version = sys.argv[1]
predict_df = pd.read_csv('infered_perturbedGX_ts_v{0}.csv'.format(version))
true_df = pd.read_csv('../signature_test_cell_4.csv')
train_data = pd.read_csv('../signature_train_cell_1_v{0}.csv'.format(version))
dev_data = pd.read_csv('../signature_dev_cell_1.csv')
test_data = pd.read_csv('../signature_test_cell_1.csv')
exclude_cells = set(dev_data.cell_id) | set(test_data.cell_id)

assert predict_df.shape == true_df.shape, "two dataset has different shape"

corre_result = []

for i in range(len(predict_df)):
    corre_result.append(pearsonr(predict_df.iloc[i, 5:], true_df.iloc[i, 5:])[0])

new_data = predict_df[np.array(corre_result) > 0.5]

concat_data = pd.concat([train_data, new_data])
concat_data = concat_data[concat_data['sig_id'].str.contains('24H') & concat_data['pert_idose'].isin(["0.04 um", "0.12 um", "0.37 um", "1.11 um", "3.33 um", "10.0 um"])]
unique_data = concat_data.drop_duplicates(new_data.columns[:5])
# unique_data = unique_data[~unique_data['cell_id'].isin(exclude_cells)]
new_version = int(version) + 1

print(len(train_data), len(unique_data))
unique_data.to_csv('../signature_train_cell_1_v{0}.csv'.format(new_version), index = False)
