# -*- coding: utf-8 -*-
"""
@author: liamo
"""
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from cobra.io import read_sbml_model

ROOT_FOLDER = "C:\\Users\\liamo\\Documents\\BIOINF\\PRETHESIS\\DATA\\"
metadata_clean_path = os.path.join(ROOT_FOLDER,"metadata_clean.csv")
# GTEX_META_PATH = os.path.join(ROOT_FOLDER,"GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
# tissue_FILTER_PATH = os.path.join(ROOT_FOLDER,"gtexv8_4tissues_rin6_tissue_list_v2.csv")
MODEL_PATH = os.path.join("C:\\Users\\liamo\\Documents\\BIOINF\\PRETHESIS","Human-GEM","model", "Human-GEM_ct_feb21.xml")
GTEX_CLEAN10_PATH = os.path.join(ROOT_FOLDER,"gtexv8_4tissues_rin6_tpm_log10_localT2.csv")

# import matlab.engine
# h1_tinit_path = os.path.join(ROOT_FOLDER,"Human1_Publication_Data_Scripts","Human1_Publication_Data_Scripts",
#                              "tINIT_GEMs","run_tINIT_outputs","GTEx","tINIT_GTEx_outputs.mat")
# eng = matlab.engine.start_matlab()
# h1_tinit = eng.load(h1_tinit_path) #too big do not do this
# eng.quit()

tissue_meta = pd.read_csv(metadata_clean_path,index_col=0,header=0)
data_clean = pd.read_csv(GTEX_CLEAN10_PATH,header=0,index_col=[0,1])
# ensids = {k:k.split(".")[0] for k in list(data_clean.columns)}

###
# print(tissue_meta["AGE"].groupby(tissue_meta["Tissue"]).value_counts().unstack(),end="\n\n")
#metadata table
pd.concat((tissue_meta["AGE"].groupby(tissue_meta["Tissue"]).value_counts().unstack(),
           tissue_meta["SEX"].groupby(tissue_meta["Tissue"]).value_counts().unstack())
           ,axis=1,keys=["AGE","SEX"])
# 3 age groups of 20 ys
pd.concat((tissue_meta["AGE_NEW"].groupby(tissue_meta["Tissue"]).value_counts().unstack(),
           tissue_meta["SEX"].groupby(tissue_meta["Tissue"]).value_counts().unstack())
           ,axis=1,keys=["AGE","SEX"])

# diff = set(meta_data["SUBJID"]).difference(set(donors))
# meta_data = meta_data.loc[~meta_data["SUBJID"].isin(diff),:]
# print(meta_data["SEX"].value_counts())
# print(meta_data["AGE"].value_counts())

# data_clean.index = pd.MultiIndex.from_tuples(
#     list(zip(tissues_clean,data_clean.index)), names=('tissue', 'sample'))

# for t in tissues_clean.unique():
#     print(data_clean.iloc[
#           tissues_clean[tissues_clean==t].index,:])

#only expressed
data_exp = data_clean.loc[:,data_clean.mean()>=5*np.log10(2)]

from sklearn.decomposition import PCA
pca_data = PCA(n_components=0.95)
pca_data.fit(data_exp)
x_reduced = pca_data.transform(data_exp)

i95 = np.where(pca_data.explained_variance_ratio_.cumsum()*100 >= 95)[0][0]
fig, ax = plt.subplots(1,2,figsize=(13, 5))
ax[0].bar(range(max(2,i95)), pca_data.explained_variance_ratio_[:max(2,i95)]*100)
ax[0].set_xticks(range(max(2,i95)))
ax[0].set_xticklabels(['PC'+str(i) for i in range(1,max(3,i95+1))], rotation=90,fontsize=8)
ax[0].set_title("Explained Variance by PC")
ax[0].set_ylabel("Percentage")

for t in tissues_clean.unique():
    sp = tissues_clean.index[tissues_clean==t]-1
    ax[1].plot(x_reduced[sp,0],x_reduced[sp,1], 'o' , label=t)
ax[1].set_title("PCA")
ax[1].legend(title="Tissue",bbox_to_anchor=(1.05, 1), loc='upper left')
ax[1].set_xlabel('PC1')
ax[1].set_ylabel('PC2')
plt.show()

### Tissue ###
tissue_df = pd.DataFrame(np.zeros((data_clean.shape[1],3)),index=list(data_clean.columns),
                         columns=["res_Breast_vs_Adipose","res_Kidney_vs_Adipose","res_Liver_vs_Adipose"])
for x in ["res_Breast_vs_Adipose","res_Kidney_vs_Adipose","res_Liver_vs_Adipose"]:
    file = pd.read_csv(os.path.join(ROOT_FOLDER,x+".csv"),header=0,index_col=0)
    for index, row in file.iterrows(): tissue_df.loc[index,x] = row["log2FoldChange"]
tissue_df.rename(columns=dict(zip(tissue_df.columns,
                                  ["BA","KA","LA"])),inplace=True)
# tissue_df.rename(index=ensids,inplace=True)

(tissue_df != 0).sum()
#de only
tissue_df.loc[(tissue_df != 0).any(axis=1),:].index


### Age ###
age_df = pd.DataFrame(np.zeros((data_clean.shape[1],2)),index=list(data_clean.columns),
                         columns=["res_40.59_vs_20.39","res_60.79_vs_20.39"])
for x in ["res_40.59_vs_20.39","res_60.79_vs_20.39"]:
    file = pd.read_csv(os.path.join(ROOT_FOLDER,x+".csv"),header=0,index_col=0)
    for index, row in file.iterrows(): age_df.loc[index,x] = row["log2FoldChange"]
age_df.rename(columns=dict(zip(tissue_df.columns,
                                  ["40-59vs20-39","60-79vs20-39"])),inplace=True)
# age_df.rename(index=ensids,inplace=True)

### Gender ###
#1=M, 2=F
gender_df = pd.read_csv(os.path.join(ROOT_FOLDER,"res_SEX_2_vs_1.csv"),header=0,index_col=0) #.sort_values(by=["padj"],axis=0)
gender_df = gender_df.drop(["padj"],axis=1)
# gender_df.rename(index=ensids,inplace=True)
gender_df.rename(columns={"log2FoldChange":"2vs1"},inplace=True)
# gender_df.index.duplicated().sum()

de_df = pd.concat((tissue_df,gender_df,age_df),axis=1,keys=["TISSUE","SEX","AGE"])
# del age_df,sex_df,tissue_df

# from cobra.core.gene import Gene, ast2str, eval_gpr, parse_gpr
model = read_sbml_model(MODEL_PATH)
gprs = {k.id:[r.id for r in k.reactions] for k in list(model.genes)}

###
#Sex-biased genes GTEx dataset
sb_path = os.path.join(ROOT_FOLDER,"GTEx_Analysis_v8_sbgenes","signif.sbgenes.txt")
sb_genes = pd.read_csv(sb_path,sep="\t",header=0) #104864 x 5
#filter to tissues of interest
sb_genes = sb_genes.loc[sb_genes["tissue"].isin(
    ["Breast_Mammary_Tissue","Kidney_Cortex","Liver","Adipose_Subcutaneous"]),:] #6351
#significant LFSR (≤ 0.05) -  already filtered
sb_genes["effsize"].describe()
# sb gene count per tissue
sb_genes.groupby(sb_genes["tissue"]).count()["gene"]

len(set(data_clean.columns).intersection(set(sb_genes["gene"]))) #3871

# tissue_active = {}
# for k in tissue_ids.keys():
#     tgenes = data_exp.iloc[tissue_ids[k],:]
#     tgactive = tgenes.mean()>=5
#     tissue_active[k] = tgactive
# active = pd.DataFrame()
# active = active.from_dict(tissue_active)
# active = active[np.sum(active,axis=1)>0] #active in at least 1 tissue, 410
# del tissue_active
#
# # active.loc[np.sum(active,axis=1)==4,:].index #active in all, 117
# active.loc[np.sum(active,axis=1)==4,:].shape[0]
# # active.loc[np.sum(active,axis=1)==1,:].index #unique active, 242
# active.loc[np.sum(active,axis=1)==1,:].sum() #no unique Adipose genes
#
# print(meta_data.loc[meta_data["SEX"]==1,"AGE"].describe()) #male age
# print(meta_data.loc[meta_data["SEX"]==2,"AGE"].describe()) #female age
# #convert age into 1-6 class?
# g = [meta_data["SEX"].iloc[np.where(meta_data["SUBJID"]==x)[0][0]] for x in samples]
# gender = pd.DataFrame(list(zip(tissues["SMTSD"],g)),columns=["SMTSD","SEX"])
# gender.value_counts()
#
# from Bio import Entrez
# Entrez.email = "pg40961@alunos.uminho.pt" # Always tell NCBI who you are
# id_list = list(active.index)
# search_results = Entrez.read(Entrez.epost(db="gene", id=",".join(id_list),idtype="acc")) #or nuccore?
# webenv = search_results["WebEnv"]
# query_key = search_results["QueryKey"]
# handle = Entrez.esummary(db="gene",id=",".join(id_list),idtype="acc",
#                          webenv=webenv,query_key=query_key)
# record = Entrez.read(handle)
# # record['DocumentSummarySet']['DocumentSummary'][0].keys()
# active_sum = pd.DataFrame()
# active_sum = active_sum.from_dict(record['DocumentSummarySet']['DocumentSummary'])
# active_summary_path = os.path.join(ROOT_FOLDER,"active_gene_esummary.csv")
# active_sum.to_csv(active_summary_path,index=False)


