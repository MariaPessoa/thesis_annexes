# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 16:17:43 2020

@author: liamo
"""
import pandas as pd
import numpy as np
import os,sys
from thresholding import local2_thresholding
from sklearn.feature_selection import VarianceThreshold

#RNA-Seq preprocessing for integration
ROOT_FOLDER = "C:\\Users\\liamo\\Documents\\BIOINF\\PRETHESIS" #"/home/mpessoa/LINUX"

#more tissue
# data_read = pd.read_csv(os.path.join(ROOT_FOLDER,"DATA","gtexv8_rin6_more.csv"),header=0,index_col=0).transpose()
# new_tissue=["Adipose - Subcutaneous","Breast - Mammary Tissue","Kidney - Cortex","Liver",
#             "Whole Blood","Stomach","Lung","Brain - Cortex","Muscle - Skeletal",
#              "Pancreas","Colon - Transverse"]
             
#median TPM
# data_read = pd.read_csv(os.path.join(ROOT_FOLDER,"DATA","GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"),
# header=0,skiprows=[0,1],sep="\t",index_col="Name").drop(["Description"],axis=1).transpose()
# data_read=data_read.loc[new_tissue,:].sort_index()

# TPM
data_read = pd.read_csv(os.path.join(ROOT_FOLDER,"DATA","gtexv8_4tissues_rin6_tpm.csv"),index_col=0,header=0).transpose()

samples = pd.read_csv(os.path.join(ROOT_FOLDER,"DATA","results",
"gtl50_gtu90_lt50_ft_4tissues_rin6_prmat_csm_nodrugex.csv"),header=0,index_col=[0,1,2])
data_read2 = data_read.loc[samples.loc(axis=0)["fastcore"].index.get_level_values(0),:]
del samples

# grouped TPM
# data_read = pd.read_csv(os.path.join(ROOT_FOLDER,"DATA","gtex_v8_median_by_TSAn.csv"),header=[0,1,2],index_col=0).transpose()
# ognames = data_read.index
# flatnames = ['_'.join(col) for col in data_read.index.values]
# data_read.index = flatnames

# print("VarianceThreshold > 0")
# constant_filter = VarianceThreshold(threshold=0)
# constant_filter.fit(data_read) #exp_data
# data_read2 = data_read.loc[:,constant_filter.get_support()]
del data_read
# data_read2 = data_read
# print("var filter shape:",data_read2.shape)

qvalues = [0.1, 0.25,0.5,0.9]
quantiles = data_read2.quantile(qvalues)
maxexp = np.log(data_read2.max().max())
global_thresholds = quantiles.T.apply(lambda x: x.mean())

gtl, gtu, lt = global_thresholds.loc[0.5], global_thresholds.loc[0.9], quantiles.loc[0.5]
scrs = {}
for sample in list(data_read2.index):
    upp_activity = local2_thresholding(data_read2.loc[sample,:], gtu, gtl, lt, maxexp)
    scrs[sample] = upp_activity

gene_scores = pd.DataFrame.from_dict(scrs).T
gene_scores.to_csv(os.path.join(ROOT_FOLDER,"DATA","gtexv8_4tissues_rin6_NOvarthresh_gtl50_gtu90_lt50.csv"))

#grouped
# gene_scores.index = ognames
# gene_scores.to_csv(os.path.join(ROOT_FOLDER,"DATA","gtexv8_median_by_TSAn_gtl50_gtu90_lt50.csv"))

