import pandas as pd
from scipy.stats import pearsonr
import os
import numpy as np
import sys

version = sys.argv[1]
predict_df = pd.read_csv('low_quality_perturbedGX_v{0}.csv'.format(version))
true_df = pd.read_csv('../signature_test_cell_4.csv')
train_data = pd.read_csv('../signature_train_cell_4_v{0}.csv'.format(version))

assert predict_df.shape == true_df.shape, "tow dataset has different shape"

corre_result = []

for i in range(len(predict_df)):
    corre_result.append(pearsonr(predict_df.iloc[i, 5:], true_df.iloc[i, 5:])[0])

new_data = true_df[np.array(corre_result) > 0.5]

concat_data = pd.concat([train_data, new_data])
unique_data = concat_data.drop_duplicates(new_data.columns[:5])
new_version = int(version) + 1

print(len(train_data), len(unique_data))
unique_data.to_csv('../signature_train_cell_4_v{0}.csv'.format(new_version), False)
