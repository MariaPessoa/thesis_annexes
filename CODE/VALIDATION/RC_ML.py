# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 16:17:43 2020

@author: liamo
"""
import pandas as pd
import numpy as np
import os,sys
from ML_func import prepare_targets,feature_select,ClassificationCV #,minmax_norm
from sklearn.feature_selection import VarianceThreshold

ROOT_FOLDER = "/home/mpessoa/LINUX"
# tissue_meta = pd.read_csv(os.path.join(ROOT_FOLDER,"DATA","metadata_clean.csv"),header=0,index_col=0)
tissue_meta = pd.read_csv(os.path.join(ROOT_FOLDER,"DATA","allrin6_metadata_clean.csv"),header=0,index_col=0)

#reaction content
ft_prmat = pd.read_csv(os.path.join(ROOT_FOLDER,"results","CSMS","gtl50_gtu90_lt50_ft_4tissues_rin6_prmat_csm_nodrugex.csv"),
                          header=0,index_col=1).drop(["Unnamed: 0", "Unnamed: 2"], axis=1)
tinit_prmat = pd.read_csv(os.path.join(ROOT_FOLDER,"results","CSMS","gtl50_gtu90_lt50_tinit_4tissues_rin6_prmat_csm_nodrugex.csv"),
                          header=0,index_col=1).drop(["Unnamed: 0", "Unnamed: 2"], axis=1)
                          
ft_more = pd.read_csv(os.path.join(ROOT_FOLDER,"results","CSMS","gtl50_gtu90_lt50_ft53_more_rin6_csm_nodrugex.csv"),
                          header=0,index_col=1).drop(["Unnamed: 0", "Unnamed: 2"], axis=1)
tinit_more = pd.read_csv(os.path.join(ROOT_FOLDER,"results","CSMS","gtl50_gtu90_lt50_tinit20_more_rin6_csm_nodrugex.csv"),
                          header=0,index_col=1).drop(["Unnamed: 0", "Unnamed: 2"], axis=1)

ft_all = pd.concat((ft_prmat,ft_more))
tinit_all = pd.concat((tinit_prmat,tinit_more))

print("VarianceThreshold > 0")
constant_filter = VarianceThreshold(threshold=0)
constant_filter.fit(tinit_all) #tinit_prmat
tinit_prmat_var = tinit_all.loc[:,constant_filter.get_support()]
print("tinit filter df:",tinit_prmat_var.shape)
del tinit_prmat

constant_filter = VarianceThreshold(threshold=0)
constant_filter.fit(ft_all) #ft_prmat
ft_prmat_var = ft_all.loc[:,constant_filter.get_support()]
print("ft filter df:",ft_prmat_var.shape)
del ft_prmat

# samples = tinit_prmat_var.index
Ytin = prepare_targets(tissue_meta.loc[tinit_prmat_var.index, "Tissue"])
Yft = prepare_targets(tissue_meta.loc[ft_prmat_var.index, "Tissue"])
# Y = tissue_meta.loc[samples,"Tissue"]+tissue_meta.loc[samples,"AGE_NEW"]
# Y_en = prepare_targets(Y)


print("Reaction content ML models...")
print("tinit...")
# rc_vr = pd.DataFrame(VarianceThreshold(threshold=0).fit_transform(ft_prmat2))
exp_data_500 = pd.DataFrame(feature_select(500, tinit_prmat_var, Ytin), index=tinit_prmat_var.index) #Y_en
rc_res = pd.DataFrame(ClassificationCV(exp_data_500, Ytin, rep=20))
rc_res.to_csv(os.path.join(ROOT_FOLDER,"results","g50-90_lt50_alltinit_tissue_rc_ml_nodrugex.csv"))

print("fastcore...")
exp_data_500 = pd.DataFrame(feature_select(500, ft_prmat_var, Yft), index=ft_prmat_var.index)
rc_res = pd.DataFrame(ClassificationCV(exp_data_500, Yft, rep=20))
rc_res.to_csv(os.path.join(ROOT_FOLDER,"results","g50-90_lt50_allft_tissue_rc_ml_nodrugex.csv"))
