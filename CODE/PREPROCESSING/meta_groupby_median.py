# -*- coding: utf-8 -*-
"""
@author: liamo
"""
import os,sys
import pandas as pd
import numpy as np

ROOT_FOLDER = "/home/mpessoa/LINUX"
# ROOT_FOLDER = "C:\\Users\\liamo\\Documents\\BIOINF\\PRETHESIS"

#TPM
GTEX_DATA_PATH = os.path.join(ROOT_FOLDER,"DATA","GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct")
#metadata
GTEX_ATT_PATH = os.path.join(ROOT_FOLDER,"DATA","GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
GTEX_META_PATH = os.path.join(ROOT_FOLDER,"DATA","GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

# data_atts = pd.read_csv(GTEX_ATT_PATH,header=0,sep="\t",usecols=["SAMPID","SMTSD","SMRIN"])
# meta_data = pd.read_csv(GTEX_META_PATH,header=0,sep="\t")
data_read = pd.read_csv(GTEX_DATA_PATH,header=2,sep="\t",index_col="Name").drop(["Description"],axis=1)

# donors = ["-".join(x.split("-")[:2]) for x in data_atts.loc[:,"SAMPID"]]
# tissue_meta = pd.DataFrame(list(data_atts["SMTSD"]),
                           # index=list(data_atts["SAMPID"]),columns=["Tissue"])

# donors_indx = [np.where(meta_data["SUBJID"]==x)[0][0] for x in donors]
# tissue_meta["SEX"] = [meta_data["SEX"].iloc[x] for x in donors_indx]
# tissue_meta["AGE"] = [meta_data["AGE"].iloc[x] for x in donors_indx]

#6 age groups -> 3 groups
# new_age_gr = {"20-29":"20-39","30-39":"20-39",
              # "40-49":"40-59","50-59":"40-59",
              # "60-69":"60-79","70-79":"60-79"}
# tissue_meta["AGE_NEW"] = [new_age_gr[x] for x in list(tissue_meta["AGE"])]
# tissue_meta.to_csv(os.path.join(ROOT_FOLDER,"DATA","all_metadata_clean.csv")) #gtexv8_metadata_clean

tissue_meta = pd.read_csv(os.path.join(ROOT_FOLDER,"DATA","all_metadata_clean.csv"),header=0,index_col=0)

#no filter
supersamp = data_read.groupby([tissue_meta["Tissue"],tissue_meta["SEX"],tissue_meta["AGE_NEW"]],axis=1).median()
print(supersamp.head(2))
supersamp.to_csv(os.path.join(ROOT_FOLDER,"DATA","gtexv8_median_by_TSAn.csv"))

