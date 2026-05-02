#!/bin/bash

python collect_low_conf_data.py \
--file_with_confidence_ph1 /workspace/l1000_ph2/Bayesian_GSE92742_with_confidence.csv \
--file_with_confidence_ph2 /workspace/l1000_ph2/Bayesian_GSE70138_with_confidence.csv \
--side_effect_data /workspace/DeepCE/DeepCE/data/side_effect/SIDER_PTs.csv \
--strongest_sig_file /workspace/DeepCE/DeepCE/data/side_effect/meta_LINCS_strongest_signatures.csv \
--phase_data_file_ph1 /workspace/l1000_ph2/Bayesian_GSE92742_Level5_COMPZ_n361481x978.h5 \
--phase_data_file_ph2 /workspace/l1000_ph2/Bayesian_GSE70138_Level5_COMPZ_n116218x978.h5 \
--drugs_filter_file /workspace/DeepCE/DeepCE/data/drug_smiles_new.csv \
--CellLineDGX_input_file_persistance /workspace/DeepCE/DeepCE/data/side_effect/SIDER_PTs_0.3_CellLineDGX.csv
