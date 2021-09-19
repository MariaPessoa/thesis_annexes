# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 16:17:43 2020

@author: liamo
"""
import pandas as pd
import os,sys
from cobra.io import read_sbml_model
from cobra.flux_analysis import find_blocked_reactions

ROOT_FOLDER = "/home/mpessoa/LINUX"

print("Reading model...")
model = read_sbml_model(os.path.join(ROOT_FOLDER, 'Human-GEM_ct_nodrugex.xml'))
print("Reading data...")
csms = pd.read_csv(os.path.join(ROOT_FOLDER,"results","CSMS",
                          "gtl50_gtu90_lt50_ft_4tissue_novar_1tpmthresh_csm_nodrugex.csv"),
                          header=0,index_col=[1]).drop(["Unnamed: 0","Unnamed: 2"], axis=1)
print("csms shape:",csms.shape)

with model as context_specific_model:
    a = 1
    for i,row in csms.iterrows():
        print(a,"out of",csms.shape[0])
        to_remove = row[(row==False)].index.to_list()
        for rid in to_remove: context_specific_model.reactions.get_by_id(rid).knock_out()
        # context_specific_model.remove_reactions(to_remove)

        blocked = find_blocked_reactions(context_specific_model)
        # print(i," - blocked:",len(blocked))
        total = csms.loc[i,:].sum()
        csms.loc[i,blocked] = False
        # print("Before:",total,"After:",csms.loc[i,:].sum())
        a+=1

csms.to_csv(os.path.join(ROOT_FOLDER,"results","CSMS",
                          "gtl50_gtu90_lt50_ft_4tissue_novar_1tpmthresh_csm_noblock_nodrugex.csv"))

