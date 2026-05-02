import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error
import numpy as np

def get_pearsonr_var_distance(input_series):
    truth_series = input_series['ehill_truth']
    pred_series = input_series['ehill_pred']
    if len(truth_series) < 6:
        return pd.Series({'pearson': None, 'spearman': None, 'variance': None, 'mse': None})
    return pd.Series({'pearson': pearsonr(truth_series.values.reshape(-1), pred_series.values.reshape(-1))[0], 
            'spearman': spearmanr(truth_series.values.reshape(-1), pred_series.values.reshape(-1))[0], 
            'variance': input_series['ehill_truth'].var(), 
            'mse': mean_squared_error(truth_series.values.reshape(-1), pred_series.values.reshape(-1))})

if __name__ == "__main__":
    
    pred_ehill = pd.read_csv('predicted_ehill_results.csv')
    truth_ehill = pd.read_csv('all_cleaned_ehill_data_test.csv')
    merged_ehill = pred_ehill.merge(truth_ehill, on = ['sig_id', 'pert_id', 'pert_type', 'cell_id', 'pert_idose'], suffixes = ('_pred', '_truth'))

    result_df = merged_ehill.groupby(['pert_id', 'pert_type', 'cell_id']).apply(get_pearsonr_var_distance).reset_index()

    drug_wise_pearson_mean_var = result_df.groupby('pert_id')[['pearson','spearman']].agg(['mean', 'var']).reset_index()
    cell_wise_pearson_mean_var = result_df.groupby('cell_id')[['pearson','spearman']].agg(['mean', 'var']).reset_index()

    drug_wise_pearson_mean_var.to_csv('drugwise_stats.csv', index = False)
    cell_wise_pearson_mean_var.to_csv('cellwise_stats.csv', index = False)

    var_threshold, mse_threshold = np.quantile(result_df['variance'].dropna(), 0.75), np.quantile(result_df['mse'].dropna(), 0.1)
    sel_drug_cell = result_df[(result_df['variance'] > var_threshold) & (result_df['mse'] < mse_threshold)]

    sel_ehill_data = merged_ehill.merge(sel_drug_cell[['pert_id', 'pert_type', 'cell_id']], on = ['pert_id', 'pert_type', 'cell_id'])
    sel_ehill_data.to_csv('sel_ehill_data.csv', index=False)