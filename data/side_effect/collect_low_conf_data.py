'''
This script is used to extract a signature list of sig_id which have low confidence and also strongest signals
'''

import pandas as pd
import argparse
from L1000_data_builder import L1000DataBuilder, ConfidenceDataBuilder
import pdb

dosage_map = {'1.11 µM':'1.11 um', 
    '40 µM': '40 um', '0.04 um': '0.04 um', '0.37 um': '0.37 um', 
    '100 nM': '0.1 um', '1 µM': '1 um', '0.04 µM': '0.04 um', '10.05 um': '10.0 um', 
    '1 nM': '1 nm', '20 µM': '20 um', '3.33 µM': '3.33 um', '1.11 um': '1.11 um', 
    '5 µM': '5 um', '3.35 um': '3.33 um', '3 µM': '3 um', '100 µM': '100 um', 
    '10 µM': '10.0 um', '50 µM': '50 um', '80 µM': '80 um', '30 µM': '30 um', 
    '90 µM': '90 um', '0.12 um': '0.12 um', '500 nM': '0.5 um', '1.12 um': '1.11um', 
    '10 nM': '10.0 nm', '10.0 um': '10.0 um', '20.0 um': '20.0 um', 
    '3.33 um': '3.33 um', '0.12 µM': '0.12 um', '0.37 µM': '0.37 um'}

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    ## use Bayesian_***_with_confidence.csv
    parser.add_argument("--file_with_confidence_ph1", help="the file directory of ph1 data with confidence")
    parser.add_argument("--file_with_confidence_ph2", help="the file directory of ph2 data with confidence")
    ## use either FAERS_offsides_PTs.csv or SIDER_PTs.csv
    parser.add_argument("--side_effect_data", help="the file directory of side effect data")
    ## use meta_LINCS_strongest_signatures.csv
    parser.add_argument("--strongest_sig_file", help="the file to save pert_id map to the corresponding strongest signature")
    ## use Bayesian_GSE***_Level5_COMPZ_n116218x978.h5
    parser.add_argument("--phase_data_file_ph1", help = 'the file to save the real phase 1 data (the signatures)')
    parser.add_argument("--phase_data_file_ph2", help = 'the file to save the real phase 2 data (the signatures)')
    parser.add_argument("--drugs_filter_file", help = 'the file to help to select drugs')
    ## the file directory to save the CellLineDGX/DeepCE model input file
    parser.add_argument("--CellLineDGX_input_file_persistance", help = "the file directory to save the CellLineDGX/DeepCE model input file")

    args = parser.parse_args()
    confidence_file_ph1 = args.file_with_confidence_ph1
    confidence_file_ph2 = args.file_with_confidence_ph2
    side_effect_file = args.side_effect_data
    strongest_sig_file = args.strongest_sig_file
    phase_data_file_ph1 = args.phase_data_file_ph1
    phase_data_file_ph2 = args.phase_data_file_ph2
    drugs_filter_file = args.drugs_filter_file
    CellLineDGX_input_file_persistance = args.CellLineDGX_input_file_persistance

    confidence_data_builder = ConfidenceDataBuilder(confidence_file_ph1, confidence_file_ph2, datatype = 'file')
    phase_sig_data_builder = L1000DataBuilder(phase_data_file_ph1, phase_data_file_ph2, phase_data_type= 'h5_file')

    drugs_names = set(pd.read_csv(drugs_filter_file).broad_cpd_id)
    confidence_df_ph1 = confidence_data_builder.retrive_confidence_data(columns=['sig_id','pert_id','cell_id','pert_idose','pert_type','confidence'], phase_num = 1)
    confidence_df_ph2 = confidence_data_builder.retrive_confidence_data(columns=['sig_id','pert_id','cell_id','pert_idose','pert_type','confidence'], phase_num = 2)
    confidence_df_ph1 = confidence_df_ph1[confidence_df_ph1.cell_id.isin(['A549', 'MCF7', 'HCC515', 
                        'HEPG2', 'HS578T', 'PC3', 'SKBR3', 'MDAMB231', 'JURKAT', 'A375', 'BT20', 'HELA', 'HT29', 'HA1E', 'YAPC'])]
    confidence_df_ph1.pert_idose = confidence_df_ph1.pert_idose.apply(lambda x: dosage_map[x] if x in dosage_map else None)
    confidence_df_ph1 = confidence_df_ph1[confidence_df_ph1.pert_idose.isin(["0.04 um", "0.12 um", "0.37 um", "1.11 um", "3.33 um", "10.0 um"])]
    confidence_df_ph1 = confidence_df_ph1[confidence_df_ph1.pert_id.isin(drugs_names)]
    confidence_df_ph2 = confidence_df_ph2[confidence_df_ph2.cell_id.isin(['A549', 'MCF7', 'HCC515', 
                        'HEPG2', 'HS578T', 'PC3', 'SKBR3', 'MDAMB231', 'JURKAT', 'A375', 'BT20', 'HELA', 'HT29', 'HA1E', 'YAPC'])]
    confidence_df_ph2.pert_idose = confidence_df_ph2.pert_idose.apply(lambda x: dosage_map[x] if x in dosage_map else None)
    confidence_df_ph2 = confidence_df_ph2[confidence_df_ph2.pert_idose.isin(["0.04 um", "0.12 um", "0.37 um", "1.11 um", "3.33 um", "10.0 um"])]
    confidence_df_ph2 = confidence_df_ph2[confidence_df_ph2.pert_id.isin(drugs_names)]
    ### only keep the lower confidence data when there are repeat
    confidence_df_ph1 = confidence_df_ph1.sort_values(['pert_id','cell_id','pert_idose','pert_type',
                                        'confidence'], ascending=True).groupby(['pert_id','cell_id','pert_idose','pert_type']).head(1)
    confidence_df_ph2 = confidence_df_ph2.sort_values(['pert_id','cell_id','pert_idose','pert_type',
                                        'confidence'], ascending=True).groupby(['pert_id','cell_id','pert_idose','pert_type']).head(1)
    print("confidence data are ready")
    side_effect_df = pd.read_csv(side_effect_file)
    ### pert_id,ctrl distil_ids,pert distil_ids,strongest sig_id 
    strongest_sig_df = pd.read_csv(strongest_sig_file)[['pert_id', 'strongest sig_id']]
    print("strongest sig data are ready")
    phase_data_df_ph1 = phase_sig_data_builder.retrive_phase_data(columns = None, phase_num=1)
    phase_data_df_ph2 = phase_sig_data_builder.retrive_phase_data(columns = None, phase_num=2)
    print("phase sig data are ready")

    #### 1. Find the sigid with pertid in side effect df and have low confidence
    side_confidence_ph1 = confidence_df_ph1.merge(side_effect_df, left_on = ['pert_id'], right_on = ['pert_id'])
    side_confidence_ph2 = confidence_df_ph2.merge(side_effect_df, left_on = ['pert_id'], right_on = ['pert_id'])
    side_confidence = pd.concat([side_confidence_ph1, side_confidence_ph2])
    print("merge data are ready")
    side_low_confidence = side_confidence[side_confidence['confidence'] < 0.3]
    side_low_confidence_time = side_low_confidence[side_low_confidence['sig_id'].apply(lambda x: '24H' in x)].drop('confidence', axis = 1)

    ### 2. build the input dataset for CellLineDGX model
    ### sig_id (24H),pert_id,pert_type,cell_id,pert_idose,DDR1, etc...
    gene_names = list(phase_data_df_ph1.columns)
    CellLineDGX_input_ph1 = phase_data_df_ph1.merge(side_low_confidence_time[['sig_id','pert_id','cell_id','pert_idose','pert_type']], left_index = True, right_on = 'sig_id')
    CellLineDGX_input_ph2 = phase_data_df_ph2.merge(side_low_confidence_time[['sig_id','pert_id','cell_id','pert_idose','pert_type']], left_index = True, right_on = 'sig_id')
    CellLineDGX = pd.concat([CellLineDGX_input_ph1, CellLineDGX_input_ph2])[['sig_id','pert_id','pert_type','cell_id','pert_idose'] + gene_names]
    CellLineDGX = CellLineDGX.drop_duplicates(['pert_id','cell_id','pert_idose','pert_type'])
    ### 3. persist the CellLineDGX inputfile
    print(len(set(CellLineDGX.pert_id)))
    CellLineDGX.to_csv(CellLineDGX_input_file_persistance, index = False)
    







    
