import pandas as pd
import argparse
from scipy.stats import pearsonr, spearmanr

parser = argparse.ArgumentParser(description='average_pearson')
parser.add_argument('--input_df')
args = parser.parse_args()
input_df = args.input_df

df = pd.read_csv(input_df, index_col = 0)
corr_scores = []
def get_average_pearson(df):
    for idx_i in df.index:
        for idx_j in df.index:
            if idx_i != idx_j:
                corr_scores.append(pearsonr(df.loc[idx_i,:], df.loc[idx_j,:])[0])
print(sum(corr_scores)/len(corr_scores))
