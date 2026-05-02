'''
split randomly to train, dev and test dataset
'''

import pandas as pd
import argparse
import random

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='DeepCE PreTraining')
    parser.add_argument('--file_name')
    args = parser.parse_args()
    file_name = args.file_name
    df = pd.read_csv(file_name)
    idx = list(df.index)
    random.seed(42)
    random.shuffle(idx)
    train_idx, test_idx = idx[len(idx)//5:], idx[:len(idx)//5]
    train_idx, dev_idx = train_idx[len(train_idx)//5:], train_idx[:len(train_idx)//5]
    file_base = file_name.rsplit('.cs')[0]
    ### make sure the first column has '24H' key word
    df.iloc[:, 0] = df.iloc[:, 0] + '_24H'
    df['pert_type'] =  'trt_cp'
    df = df[['sig_id', 'pert_id', 'pert_type', 'cell_id', 'pert_idose', 'ehill']]
    df.loc[train_idx, :].to_csv(file_base + '_train.csv', index = False)
    df.loc[dev_idx, :].to_csv(file_base + '_dev.csv', index = False)
    df.loc[test_idx, :].to_csv(file_base + '_test.csv', index = False)
