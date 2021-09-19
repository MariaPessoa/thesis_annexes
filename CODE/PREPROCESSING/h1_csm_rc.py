# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 16:17:43 2020

@author: liamo
"""
import pandas as pd
import numpy as np
import os,sys
import matplotlib.pyplot as plt
from scipy.io import loadmat
from cobra.io.mat import from_mat_struct
from cobra.io import read_sbml_model

ROOT_FOLDER = "C:\\Users\\liamo\\Documents\\BIOINF\\PRETHESIS"

#MATLAB to reaction content
h1_csms = loadmat(os.path.join(ROOT_FOLDER,"DATA","Human1_Publication_Data_Scripts","tINIT_GEMs","run_tINIT_outputs",
                               "GTEx","tINIT_GTEx_outputs.mat"))["INIT_output"]

name_ids = pd.DataFrame([x[0][0] for x in h1_csms[0, 0][0]])
# new_tissues = ["pancreas","stomach","lung","brain","muscle","blood","colon"]
# all_models = []
# for t in sorted(new_tissues):
#     model = from_mat_struct(h1_csms[0, 0][1][np.where(name_ids == t)[0][0]][0], model_id=t)
#     all_models.append(model)

ogmodel = read_sbml_model(os.path.join(ROOT_FOLDER,"Human-GEM","model","Human-GEM.xml"))
all_reactions = [x.id for x in ogmodel.reactions]

# h1_csm_rc = pd.DataFrame(np.zeros((len(all_reactions), len(new_tissues))),index=all_reactions,
#                           columns=sorted(new_tissues))
h1_csm_rc = pd.DataFrame(np.zeros((len(all_reactions), len(name_ids))),index=all_reactions,
                          columns=sorted(name_ids))

for i,t in enumerate(sorted(name_ids)): #sorted(new_tissues)
    print(t,i+1,"out of",len(name_ids))
    model = from_mat_struct(h1_csms[0, 0][1][np.where(name_ids == t)[0][0]][0], model_id=t)
    rcts = [r.id for r in model.reactions if r.bounds != (0,0)] #all_models[i]
    h1_csm_rc.loc[rcts,t] = 1

h1_csm_rc = h1_csm_rc.astype(bool).transpose()
h1_csm_rc.to_csv(os.path.join(ROOT_FOLDER,"DATA","results","tINIT_GTEx_all_rc.csv"))

