# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 16:17:43 2020

@author: liamo
"""
import pandas as pd
import numpy as np
import os,sys
from itertools import product
from thresholding import local2_thresholding,global_thresholding,local2_richelle
from ML_func import prepare_targets,feature_select,ClassificationCV
from sklearn.feature_selection import VarianceThreshold

ROOT_FOLDER = "/home/mpessoa/LINUX"
#gene x samples -T-> samples x genes
#new tissues
tissue_meta = pd.read_csv(os.path.join(ROOT_FOLDER,"DATA","allrin6_metadata_clean.csv"),index_col=0,header=0)

# new_tissues = ["Whole Blood","Stomach","Lung","Brain - Cortex","Muscle - Skeletal",
             # "Pancreas","Colon - Transverse"]
# save_samples_new = tissue_meta.loc[tissue_meta["Tissue"].str.contains(
    # "|".join(new_tissues))].index.to_list()
# data_read = pd.read_csv(os.path.join(ROOT_FOLDER,"DATA","GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct"),
# header=2,sep="\t", usecols=["Name"]+save_samples_new)
# data_read.index = data_read["Name"]
# data_read = data_read.drop(["Name"],axis=1)
# data_read.to_csv(os.path.join(ROOT_FOLDER,"DATA","gtexv8_rin6_more.csv"))

data_read = pd.read_csv(os.path.join(ROOT_FOLDER,"DATA","gtexv8_rin6_more.csv"),header=0,index_col=0).transpose()

# print("VarianceThreshold > 0") #do not use
# constant_filter = VarianceThreshold(threshold=0)
# constant_filter.fit(data_read) #exp_data
# exp_data2 = data_read.loc[:,constant_filter.get_support()] #data_read
# del data_read
# print("exp_data2 shape:",exp_data2.shape)
###
exp_data2=data_read

Y = tissue_meta.loc[exp_data2.index, "Tissue"] #"Tissue"
Y_en = prepare_targets(Y)                                  

qvalues = [0.1, 0.25, 0.5, 0.75, 0.9]
quantiles = exp_data2.quantile(qvalues)
maxexp = np.log(exp_data2.max().max())

global_thresholds = quantiles.T.apply(lambda x: x.mean())
global_lt2_params = list(product(range(len(qvalues)),
                                 list(zip(*np.where(np.fromfunction(lambda i, j: i < j, [len(qvalues)] * 2))))))

print("Omics data...")
exp_data_500 = pd.DataFrame(feature_select(500, exp_data2, Y_en), index=exp_data2.index) #exp_data2, exp_data2.index
omics_results = pd.DataFrame(ClassificationCV(exp_data_500, Y_en, rep=20)) #def cv=5
omics_results.to_csv(os.path.join(ROOT_FOLDER,"results","threshold_test_omics_alltissue_rin6.csv"))

print("Gene scores...")
print("\tGlobal thresholds...")
### global
results_global = []
for q in qvalues:
    current = {}
    for sample in list(exp_data2.index): #list(exp_data2.index)
        upp_activity = global_thresholding(exp_data2.loc[sample, :], global_thresholds.loc[q], 0, 0, maxexp)
        current[sample] = upp_activity

    global_df = pd.DataFrame.from_dict(current).T
    exp_data_500 = pd.DataFrame(feature_select(500, global_df, Y_en), index=global_df.index)
    results_global.append(pd.DataFrame(ClassificationCV(exp_data_500, Y_en, rep=20)))  #cv x rep

global_res = pd.concat([x.mean(axis=1) for x in results_global], axis=1,keys=qvalues)
global_res.to_csv(os.path.join(ROOT_FOLDER,"results","threshold_test_global_alltissue_rin6.csv"))

print("\tLocal T2 thresholds...")
### local t2
results_local = []
for v, k in global_lt2_params:
    current = {}
    gtl, gtu, lt = global_thresholds.iloc[k[0]], global_thresholds.iloc[k[1]], quantiles.iloc[v, :]
    for sample in list(exp_data2.index): #list(exp_data2.index)
        upp_activity = local2_thresholding(exp_data2.loc[sample,:], gtu, gtl, lt, maxexp)
        current[sample] = upp_activity

    local_df = pd.DataFrame.from_dict(current).T
    exp_data_500 = pd.DataFrame(feature_select(500, local_df, Y_en), index=local_df.index)
    results_local.append(pd.DataFrame(ClassificationCV(exp_data_500, Y_en, rep=20)))  # cv x rep

global_lt2_params2 = [(qvalues[x[0]],(qvalues[x[1][0]],qvalues[x[1][1]])) for x in global_lt2_params]
local_res = pd.concat([x.mean(axis=1) for x in results_local], axis=1,keys=global_lt2_params2)
local_res.to_csv(os.path.join(ROOT_FOLDER,"results","threshold_test_local_alltissue_rin6.csv"))

