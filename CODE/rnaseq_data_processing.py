# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 16:17:43 2020

@author: liamo
"""
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

#RNA-Seq preprocessing for integration
ROOT_FOLDER = "C:\\Users\\liamo\\Documents\\BIOINF\\PRETHESIS"

#TPM for model integration, unnormalized reads for differential expression analysis
GTEX_DATA_RAW_PATH = os.path.join(ROOT_FOLDER,"DATA","GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
GTEX_DATA_PATH = os.path.join(ROOT_FOLDER,"DATA","GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct")
#metadata for filtering
GTEX_ATT_PATH = os.path.join(ROOT_FOLDER,"DATA","GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
GTEX_META_PATH = os.path.join(ROOT_FOLDER,"DATA","GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

#output paths
GTEX_FILTER_RAW_PATH = os.path.join(ROOT_FOLDER,"DATA","gtexv8_4tissues_rin6_reads_ENS18k.csv")
# GTEX_FILTER_PATH = os.path.join(ROOT_FOLDER,"DATA","gtexv8_4tissues_rin6_tpm.csv")
# tissue_FILTER_PATH = os.path.join(ROOT_FOLDER,"DATA","gtexv8_4tissues_rin6_tissue_list_v2.csv")
GTEX_CLEAN10_PATH = os.path.join(ROOT_FOLDER,"DATA","gtexv8_4tissues_rin6_tpm_log10_localT2.csv")

data_sample = pd.read_csv(GTEX_DATA_PATH,header=2,sep="\t", nrows=5)
data_atts = pd.read_csv(GTEX_ATT_PATH,header=0,sep="\t")
meta_data = pd.read_csv(GTEX_META_PATH,header=0,sep="\t")

tissues = ["Adipose - Subcutaneous","Breast","Kidney - Cortex","Liver"]
#only 4 "kidney - medulla", which will be discontinued

#tissue sample & rna quality filter  
data = data_atts[(data_atts["SMTSD"].str.contains(
    "|".join(tissues),case=False, regex=True)) & (
        data_atts["SMRIN"] >=6)].sort_values("SMRIN",axis=0,ascending=False) #(1369, 63)

#cross checking donor ids    
save_samples = list(set(data.iloc[:,0]).intersection(set(data_sample.columns.values))) #1228

data_read = pd.read_csv(GTEX_DATA_PATH,header=2,sep="\t", usecols=["Name"]+save_samples)
# del data,data_sample,data_atts
# data_read.info() #527MB shape: (56200, 1229)

data_read.index = data_read["Name"]
data_read = data_read.drop(["Name"],axis=1)
#filter TPM > 1
data_read = data_read.loc[data_read.mean(axis=1)>1,:] #(18109, 1228)

#metadata - for filtering and analysis
tissues_clean = data.loc[data["SAMPID"].str.contains(
    "|".join(list(data_read.columns)),case=False, regex=True),["SAMPID","SMTSD"]]
donors = ["-".join(x.split("-")[:2]) for x in tissues_clean.loc[:,"SAMPID"]]

tissue_meta = pd.DataFrame(list(tissues_clean["SMTSD"]),
                           index=list(tissues_clean["SAMPID"]),columns=["Tissue"])

donors_indx = [np.where(meta_data["SUBJID"]==x)[0][0] for x in donors]
tissue_meta["SEX"] = [meta_data["SEX"].iloc[x] for x in donors_indx]
tissue_meta["AGE"] = [meta_data["AGE"].iloc[x] for x in donors_indx]

new_age_gr = {"20-29":"20-39","30-39":"20-39",
              "40-49":"40-59","50-59":"40-59",
              "60-69":"60-79","70-79":"60-79"}
tissue_meta["AGE_NEW"] = [new_age_gr[x] for x in list(tissue_meta["AGE"])]

metadata_clean_path = os.path.join(ROOT_FOLDER,"metadata_clean.csv")
# tissue_meta.to_csv(metadata_clean_path,index=True)

######################################################
#PCA
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

pca_data = PCA(n_components=0.95)
x = StandardScaler().fit_transform(data_read.transpose())
pca_data.fit(x)
x_reduced = pca_data.transform(x)

# i95 = np.where(pca_data.explained_variance_ratio_.cumsum()*100 >= 95)[0][0]
# fig, ax = plt.subplots(1,2,figsize=(12, 5))
# ax[0].bar(range(max(2,i95)), pca_data.explained_variance_ratio_[:max(2,i95)]*100)
# ax[0].set_xticks(range(max(2,i95)))
# ax[0].set_xticklabels(['PC'+str(i) for i in range(1,max(3,i95+1))],rotation=90)
# ax[0].set_title("Explained Variance by PC")
# ax[0].set_ylabel("Percentage")

### POR CORES IGUAIS ###
fig, ax = plt.subplots(figsize=(7, 5))
for t in tissue_meta["Tissue"].unique():
    meta_indx = tissue_meta.iloc[np.where(tissue_meta["Tissue"] == t)[0],:].index
    sp = np.where(data_read.columns.str.contains("|".join(list(meta_indx)), case=False, regex=True))[0]
    ax.plot(x_reduced[sp,0],x_reduced[sp,1], 'o' , label=t)
ax.set_title("PCA")
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.legend(bbox_to_anchor=(1, 1), shadow=False) #bbox_to_anchor=(1.05, 1)
fig.tight_layout()
plt.show()

# GTEX_FILTER_RAW_PATH = os.path.join(ROOT_FOLDER,"DATA","gtexv8_4tissues_rin6_reads_ENS18k.csv")
#raw counts filter
raw_data_read = pd.read_csv(GTEX_DATA_RAW_PATH,header=2,sep="\t", usecols=["Name"]+save_samples)
raw_data_read.index = raw_data_read["Name"]
raw_data_read = raw_data_read.drop(["Name"],axis=1)

raw_data_read = raw_data_read.loc[raw_data_read.index.isin(data_read.index),:]
raw_data_read.to_csv(GTEX_FILTER_RAW_PATH,index=True) #genes x samples
del raw_data_read

###############################################################
### try normalization per tissue
# tissue_ids = {k:np.where(tissue_samples["SMTSD"]==k)[0] for k in tissue_samples["SMTSD"].unique()}
tissue_ids = {}
for t in tissue_meta["Tissue"].unique():
    meta_indx = tissue_meta.iloc[np.where(tissue_meta["Tissue"] == t)[0], :].index
    sp = np.where(data_read.columns.str.contains("|".join(list(meta_indx)), case=False, regex=True))[0]
    tissue_ids[t] = sp

adipose = data_read.iloc[:,tissue_ids["Adipose - Subcutaneous"]]
breast = data_read.iloc[:,tissue_ids["Breast - Mammary Tissue"]]
liver = data_read.iloc[:,tissue_ids["Liver"]]
kidney = data_read.iloc[:,tissue_ids["Kidney - Cortex"]]
# del data_read

def normalize(indf):
    data_mean = indf.mean(axis=1)
    results = []

    #overall quantile
    q25 = np.quantile(indf, 0.25)
    q75 = np.quantile(indf, 0.75)
    unexp = (data_mean <= 0)
    q25norm = ~unexp & (data_mean <= q25)
    q75norm = ~unexp & (data_mean >= q75)
    meannorm = ~unexp & (data_mean < q75) & (data_mean > q25)
    # print(unexp.sum() + q75norm.sum() + q25norm.sum() + meannorm.sum() == indf.shape[0])

    for row in range(indf.shape[0]):
        if q75norm.iloc[row]:
            # indf.iloc[row, :] = 5 * np.log10(1 + (indf.iloc[row, :] / q75))
            results.append(5 * np.log10(1 + (indf.iloc[row, :] / q75)))
        elif q25norm.iloc[row]:
            # indf.iloc[row, :] = 5 * np.log10(1 + (indf.iloc[row, :] / q25))
            results.append(5 * np.log10(1 + (indf.iloc[row, :] / q25)))
        elif meannorm.iloc[row]:
            # indf.iloc[row, :] = 5 * np.log10(1 + (indf.iloc[row, :] / data_mean.iloc[row]))
            results.append(5 * np.log10(1 + (indf.iloc[row, :] / data_mean.iloc[row])))
        else:
            # indf.iloc[row, :] = 0
            results.append(0*indf.iloc[row, :])

    #quantile per gene?
    # q25 = indf.quantile(0.25, axis=1)
    # q75 = indf.quantile(0.75, axis=1)
    # unexp = (data_mean <= 0) | (q25 <= 0) | (q75 <= 0)
    # q25norm = ~unexp & (data_mean <= q25)
    # q75norm = ~unexp & (data_mean >= q75)
    # meannorm = ~unexp & (data_mean < q75) & (data_mean > q25)
    # # check
    # # print(unexp.sum() + q75norm.sum() + q25norm.sum() + meannorm.sum() == indf.shape[0])
    #
    # # log10
    # for row in range(indf.shape[0]):
    #     if q75norm.iloc[row]:
    #         indf.iloc[row, :] = 5 * np.log10(1 + (indf.iloc[row, :] / q75.iloc[row]))
    #     elif q25norm.iloc[row]:
    #         indf.iloc[row, :] = 5 * np.log10(1 + (row / q25.iloc[row]))
    #     elif meannorm.iloc[row]:
    #         indf.iloc[row, :] = 5 * np.log10(1 + (indf.iloc[row, :] / data_mean.iloc[row]))
    #     else:
    #         indf.iloc[row, :] = 0

    fastcore = pd.DataFrame(results,index=list(indf.index),columns=list(indf.columns)).transpose()
    # fastcore = indf.transpose()  # samples x genes
    return fastcore

# adipose_clean = normalize(adipose)
# breast_clean = normalize(breast)
# liver_clean = normalize(liver)
# kidney_clean = normalize(kidney)

# GTEX_CLEAN10_PATH = os.path.join(ROOT_FOLDER,"DATA","gtexv8_4tissues_rin6_tpm_log10_globalqts.csv")
# GTEX_CLEAN10_PATH = os.path.join(ROOT_FOLDER,"DATA","gtexv8_4tissues_rin6_tpm_log10_localqts.csv")
#
# data_clean_all = pd.concat([adipose_clean.transpose(),breast_clean.transpose(),
#                             liver_clean.transpose(),kidney_clean.transpose()],axis=1,
#           keys=["Adipose","Breast","Liver","Kidney"]).transpose()
# data_clean_all.to_csv(GTEX_CLEAN10_PATH,index=True)
# del adipose_clean, breast_clean, liver_clean, kidney_clean

data_cnorm = normalize(data_read)
data_commonnorm_path = os.path.join(ROOT_FOLDER,"DATA","gtexv8_4tissues_rin6_tpm_log10_globalqts_cnorm.csv")
data_cnorm.to_csv(data_commonnorm_path,index=True)

test = []
for t in tissue_meta["Tissue"].unique():
    meta_indx = tissue_meta.iloc[np.where(tissue_meta["Tissue"] == t)[0], :].index
    sp = np.where(data_read.columns.str.contains("|".join(list(meta_indx)), case=False, regex=True))[0][:53]
    test.extend(sp)

#Richelle, A., Chiang, A. W., Kuo, C. C., & Lewis, N. E. (2019).
#Increasing consensus of context-specific metabolic models by integrating 
#data-inferred cell functions. PLoS computational biology, 15(4), e1006867.

