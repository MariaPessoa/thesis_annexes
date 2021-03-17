# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 12:16:19 2020
@author: liamo
"""
import os, sys
from cobra.io import read_sbml_model,load_matlab_model #, write_sbml_model
from cobra.flux_analysis import find_blocked_reactions
from cobra.flux_analysis import find_essential_reactions,flux_variability_analysis,\
    pfba,gapfill
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join("C:\\Users\\liamo\\Documents\\BIOINF\\PRETHESIS\\troppo-dev\\src"))
sys.path.append(os.path.join("C:\\Users\\liamo\\Documents\\BIOINF\\PRETHESIS\\cobamp-dev\\src"))
from troppo.tasks.task_io import JSONTaskIO,ExcelTaskIO

ROOT_FOLDER = "C:\\Users\\liamo\\Documents\\BIOINF\\PRETHESIS\\"
MODEL_PATH = os.path.join(ROOT_FOLDER,"Human-GEM","model", "Human-GEM_ct.xml")
MAT_MODEL_PATH = os.path.join(ROOT_FOLDER,"Human-GEM","model", "Human-GEM.mat") #original / non consistent

HUMAN1_TASKS_PATH = os.path.join(ROOT_FOLDER,"Human-GEM","data","metabolicTasks", 'metabolicTasks_Essential.xlsx')
GTEX_CLEAN10_PATH = os.path.join(ROOT_FOLDER,"DATA","gtexv8_4tissues_rin6_tpm_log10_globalqts_cnorm.csv")
metadata_clean_path = os.path.join(ROOT_FOLDER,"DATA","metadata_clean.csv")

h1_tasks = ExcelTaskIO().read_task(HUMAN1_TASKS_PATH)
tissue_meta = pd.read_csv(metadata_clean_path,header=0,index_col=0)
data_clean = pd.read_csv(GTEX_CLEAN10_PATH,header=0,index_col=[0,1])

ft_path = os.path.join(ROOT_FOLDER,"DATA","results","b212_ft_log10_csm_min_glqts_cnorm.csv")
ft_tasks_path = os.path.join(ROOT_FOLDER,"DATA","results","b212_ft_log10_tasks_min_glqts_cnorm.json")
tinit_path = os.path.join(ROOT_FOLDER,"DATA","results","b212_tinit_log10_csm_min_glqts_cnorm.csv")
tinit_tasks_path = os.path.join(ROOT_FOLDER,"DATA","results","b212_tinit_log10_tasks_min_glqts_cnorm.json")

ft = pd.read_csv(ft_path,header=0,index_col=1).drop(["Unnamed: 0","Unnamed: 2"],axis=1)
ft_tasks = pd.read_json(ft_tasks_path) #model x tasks
tinit = pd.read_csv(tinit_path,header=0,index_col=1).drop(["Unnamed: 0","Unnamed: 2"],axis=1)
tinit_tasks = pd.read_json(tinit_tasks_path)

# tinit["biomass_human"].sum()
# ft["biomass_human"].sum()

# np.sum(ft,axis=1)
# print(tissues.value_counts())
# print(pd.DataFrame(ft_tissues).value_counts())

######################################################
### Reaction count vs tasks passed

def tasks_to_df(pddf):
    df = pd.DataFrame(np.zeros((pddf.shape[0],len(pddf.iloc[0,1].values()))))
    df_indx = [x[1] for x in pddf.iloc[:,0]]
    df.index = df_indx
    #ver nomes tasks vitor
    df_cols = dict(zip(df.columns,list(pddf.iloc[0,1].keys())))
    df.rename(columns=df_cols,inplace=True)

    for row in range(pddf.shape[0]):
        line = list(pddf.iloc[row,1].values())
        line2 = []
        for lst in line: line2.append(lst[0])
        df.iloc[row,:] = line2
    return df

ft_tasks_df = tasks_to_df(ft_tasks)
tinit_tasks_df = tasks_to_df(tinit_tasks)

rcts = pd.concat((ft.sum(axis=1),tinit.sum(axis=1)),axis=1,keys=["fastcore","tinit"])
tsks = pd.concat((ft_tasks_df.sum(axis=1),tinit_tasks_df.sum(axis=1)),axis=1,keys=["fastcore","tinit"])
full = pd.concat((rcts,tsks),axis=1,keys=["reactions","tasks"])

# def_colors = plt.rcParams['axes.prop_cycle'].by_key()['color'][:4]
# tissue_color = dict(zip(sorted(tissue_meta["Tissue"].unique()),def_colors))

#gde diferença. verificar se foi msm o modelo consistente usado como template para o tinit?
fig, ax = plt.subplots() #marker = "o" color=sample_color.values() #e a legenda?
full.plot(kind="scatter",x=("reactions","fastcore"),y=("tasks","fastcore"),ax=ax,label="fastcore")
full.plot(kind="scatter",x=("reactions","tinit"),y=("tasks","tinit"),ax=ax, color="#ff7f0e",label="tinit")

ax.set_xlabel("Reaction count")
ax.set_ylabel("Tasks passed")
ax.legend(title='Algorithm')
plt.show()

tsks.describe()
# excluding tinit's outlier
tsks.drop(tsks.index[np.where(tsks["tinit"]==tsks["tinit"].min())[0][0]],axis=0).describe()

rcts.describe()
abs(rcts.describe().loc["mean","fastcore"]-rcts.describe().loc["mean","tinit"])
rcts.groupby(tissue_meta["Tissue"]).mean() #numero de reaçoes muito semelhante entre tecidos...

#diferenças das reaçoes fastcore vs tinit
ft.groupby(tissue_meta["Tissue"]).sum()
ft.sum().describe() #muito enviesado, reaçoes estao em quase todos os modelos

#filtrar reaçoes em todos os modelos (para descobrir as unicas a um tecido)
ft.loc[:,ft.sum()<ft.shape[0]] #3054
#filtered group by
filt_gb = ft.loc[:,ft.sum()<ft.shape[0]].groupby(tissue_meta["Tissue"]).sum()
#filtrar reaçoes em nenhum modelo
unique_gb = filt_gb.loc[:,filt_gb.sum()!=0] #1303
(unique_gb==0).any().sum() #71
#reaçoes que nao estao em todos os tecidos
filt_unique_gb = unique_gb.loc[:,(unique_gb==0).any()]
((filt_unique_gb==0).sum()==3).sum() #29
filt_unique_gb.loc[:,((filt_unique_gb==0).sum()==3)].sum() #apenas 1 modelo/tecido...

#reaçoes presentes em pelo menos 1/2 dos modelos, por tecido
ft_rcts = ft.groupby(tissue_meta["Tissue"]).sum() > int(53/2)
tinit_rcts = tinit.groupby(tissue_meta["Tissue"]).sum() > int(53/2)
#reacçoes do tinit q o fastcore nao tem
rcts_diff_gb = tinit_rcts & ~ft_rcts
rcts_diff_gb.sum().value_counts()
#reacçoes presentes em pelo menos um dos tecidos
rcts_diff_gb_slim = rcts_diff_gb.loc[:,rcts_diff_gb.sum()>0] #1669

# model = read_sbml_model(MODEL_PATH)
# len(set(rcts_diff_slim.columns).intersection(set(model.boundary))) #no boundary reactions
#reaçoes presentes em TODOS os tecidos
protec_tinit = pd.DataFrame(rcts_diff_gb_slim.loc[:,rcts_diff_gb_slim.sum()==4].columns,columns=["ID"]) #1481
protec_tinit_path = os.path.join(ROOT_FOLDER,"DATA","task_rcts_protected_tinit.csv")
protec_tinit.to_csv(protec_tinit_path)

#best tinit models (task evaluation)
best_tinit = tinit.iloc[np.where(tsks["tinit"]==tsks["tinit"].max())[0],:]
bt_ft = ft.iloc[np.where(tsks["tinit"]==tsks["tinit"].max())[0],:]
rcts_diff = best_tinit & ~bt_ft
rcts_diff = rcts_diff.loc[:,rcts_diff.sum()>0] #2148
rcts_diff.sum().describe()

# rcts_diff.loc[:,rcts_diff.sum()>rcts_diff.sum().median()] #1001
len(set(rcts_diff_gb_slim.columns).intersection(set(
    rcts_diff.loc[:,rcts_diff.sum()>rcts_diff.sum().median()].columns))) #1000... nope
len(set(rcts_diff.columns).intersection(set(rcts_diff_gb_slim.columns))) #1609
bt_notgb = set(rcts_diff.columns).difference(set(rcts_diff_gb_slim.columns)) #539

PROTECTED_BEST_PATH = os.path.join(ROOT_FOLDER,"DATA","task_rcts_protected_best.csv")
protec_best = pd.DataFrame(bt_notgb,columns=["ID"])
protec_best.to_csv(PROTECTED_BEST_PATH)

#comparar as reaçoes preotegidas
PROTECTED_EXTRA_PATH = os.path.join(ROOT_FOLDER,"DATA","task_rcts_protected.csv")
protec = pd.read_csv(PROTECTED_EXTRA_PATH,header=0,index_col=0)
# len(set(rcts_diff_gb_slim.columns).intersection(set(protec["ID"]))) #nenhuma!!
# len(bt_notgb.intersection(set(protec["ID"]))) #nenhuma! thats good

#testar correr com protected reactions
# ft_protec = pd.read_csv(os.path.join(ROOT_FOLDER,"DATA","results","b40_ft_log10_csm_protec.csv"))
# ft_protec_tasks = tasks_to_df(pd.read_json(os.path.join(ROOT_FOLDER,"DATA","results","b40_ft_log10_tasks_protec.json")))
# ft_protec.sum(axis=1).describe()
# ft_protec_tasks.sum(axis=1).describe() #11-21, no dice

# testar correr com protected tinit reactions
# ft_protec = pd.read_csv(os.path.join(ROOT_FOLDER,"DATA","results","b40_ft_log10_csm_protectinit.csv"))
# ft_protec_tasks = tasks_to_df(pd.read_json(os.path.join(ROOT_FOLDER,"DATA","results","b40_ft_log10_tasks_protectinit.json")))
# ft_protec.sum(axis=1).describe()
# ft_protec_tasks.sum(axis=1).describe() #todas 13

#testar correr best tinit -> todas(?) 20, ja melhorou um pouco
#testar correr protected com o tinit n=40 -> doesnt seem to help

### Structural comparison ###
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

ft_sc = StandardScaler().fit_transform(ft)
tinit_sc = StandardScaler().fit_transform(tinit)

tsne_ft = TSNE()
tsne_ti = TSNE()
points_ft = tsne_ft.fit_transform(ft_sc) #.transpose() #ft
points_ti = tsne_ti.fit_transform(tinit_sc)

fig, ax = plt.subplots(1,2,figsize=(13,5))

for t in sorted(tissue_meta["Tissue"].unique()):
    meta_indx = tissue_meta.iloc[np.where(tissue_meta["Tissue"] == t)[0], :].index
    sp = np.where(ft.index.str.contains("|".join(list(meta_indx)), case=False, regex=True))[0]
    ax[0].plot(points_ft[sp, 0], points_ft[sp, 1], 'o', label=t)
    ax[1].plot(points_ti[sp, 0], points_ti[sp, 1], 'o', label=t) #???
    # plt.plot(points_ft[sp, 0], points_ft[sp, 1], 'o', label=t)

ax[0].set_title("Fastcore")
ax[0].set_xlabel('tSNE1')
ax[0].set_ylabel('tSNE2')
fig.suptitle('Structural comparison - tSNE')
ax[1].set_title("tINIT")
ax[1].set_xlabel('tSNE1')
ax[1].set_ylabel('tSNE2')
# ax[1].legend(bbox_to_anchor=(1, 1), shadow=False)
ax[1].legend(shadow=False)
plt.show()

pca_ft= PCA(n_components=0.95)
pca_ti= PCA(n_components=0.95)

pca_ft.fit(ft_sc)
x_reduced_ft = pca_ft.transform(ft_sc)
pca_ti.fit(tinit_sc)
x_reduced_ti = pca_ti.transform(tinit_sc)

fig, ax = plt.subplots(1,2,figsize=(13,5))
#by tissue
for t in sorted(tissue_meta["Tissue"].unique()):
    meta_indx = tissue_meta.iloc[np.where(tissue_meta["Tissue"] == t)[0], :].index
    sp = np.where(ft.index.str.contains("|".join(list(meta_indx)), case=False, regex=True))[0]
    # sp = np.where(ft.index.str.contains("|".join(meta_indx),case=False, regex=True))[0]
    ax[0].plot(x_reduced_ft[sp, 0], x_reduced_ft[sp, 1], 'o', label=t)
    ax[1].plot(x_reduced_ti[sp, 0], x_reduced_ti[sp, 1], 'o', label=t)
ax[0].set_title("Fastcore")
ax[0].set_xlabel('PC1')
ax[0].set_ylabel('PC2')
fig.suptitle('Structural comparison - PCA') #, fontsize=16
ax[1].set_title("tINIT")
ax[1].set_xlabel('PC1')
ax[1].set_ylabel('PC2')
ax[1].legend(bbox_to_anchor=(1, 1), shadow=False)
plt.show()

#by sex
# for t in sorted(tissue_meta["SEX"].unique()):
#     meta_indx = tissue_meta.iloc[np.where(tissue_meta["SEX"] == t)[0], :].index
#     sp = np.where(ft.index.str.contains("|".join(list(meta_indx)), case=False, regex=True))[0]
#     # sp = np.where(ft.index.str.contains("|".join(meta_indx),case=False, regex=True))[0]
#     plt.plot(x_reduced[sp, 0], x_reduced[sp, 1], 'o', label=t)
# plt.title("Structural comparison")
# plt.xlabel('PC 1')
# plt.ylabel('PC 2')
# plt.legend(bbox_to_anchor=(1, 1), shadow=False)
# plt.show()
#
# #by age
# for t in sorted(tissue_meta["AGE_NEW"].unique()):
#     meta_indx = tissue_meta.iloc[np.where(tissue_meta["AGE_NEW"] == t)[0], :].index
#     sp = np.where(ft.index.str.contains("|".join(list(meta_indx)), case=False, regex=True))[0]
#     # sp = np.where(ft.index.str.contains("|".join(meta_indx),case=False, regex=True))[0]
#     plt.plot(x_reduced[sp, 0], x_reduced[sp, 1], 'o', label=t)
# plt.title("Structural comparison")
# plt.xlabel('PC 1')
# plt.ylabel('PC 2')
# plt.legend(bbox_to_anchor=(1, 1), shadow=False)
# plt.show()

#unique reactions
#metadata influence

sum_by_tissue = []
for t in set(ft_tissues):
    sp = np.where(ft_tissues == t)[0]
    sum_by_tissue.append(ft.iloc[sp,:].sum())
sumdf = pd.concat(sum_by_tissue,axis=1,keys=set(ft_tissues)) #11895
#remove reactions not present in any model
sumdf = sumdf.loc[np.sum(sumdf,axis=1)!=0,:] #10789 #change to any/all / .sum(axis=1)
#reactions unique to a single tissue
unique = sumdf.loc[np.sum(sumdf>0,axis=1)==1,:] #140
#present in at least 2 models of each tissue
unique.loc[np.sum(unique,axis=1)>1,:] #46


### Structural similarity ###
# tissue_ids = {k:list(np.where(ft_tissues == k)[0]) for k in set(ft_tissues)}
# kidney = ft.iloc[tissue_ids["Kidney"],:]
flatten = lambda t: [item for sublist in t for item in sublist]
#sample, n=5 / tissue
tissue_ids = [list(np.where(ft_tissues == k)[0][:5]) for k in set(ft_tissues)]
ft_slim = ft.iloc[flatten(tissue_ids), :]

from sklearn.metrics import pairwise_distances
ham_dist = pairwise_distances(ft_slim,metric="hamming")
labels = ["-".join([x.split("_")[0],x.split("-")[-1]]) for x in ft_slim.index]

fig = plt.figure(figsize=(5, 5))
plt.matshow(ham_dist, fignum=fig.number)
plt.xticks(range(ham_dist.shape[1]), labels, rotation=90) #, fontsize=12
plt.yticks(range(ham_dist.shape[1]), labels) #, fontsize=12
cb = plt.colorbar()
# cb.ax.tick_params(labelsize=14)
plt.title('Structural similarity') #, fontsize=16

plt.show()

### Subsystem coverage ###
model = read_sbml_model(MODEL_PATH)
mat_model = load_matlab_model(MAT_MODEL_PATH)
reaction_ids = [x.id for x in list(model.reactions)]
# dir(list(model.reactions)[0])
#check
# len(set([x.id for x in mat_model.reactions]).intersection(set([x.id for x in model.reactions])))
# len(set([x.id for x in model.reactions]).intersection(set([x.id for x in mat_model.reactions])))

#sbml doesnt have them: ""
subsystems =  set([x.subsystem.split("'")[1] for x in list(mat_model.reactions)])
#mesma ordem, alfabetica
model_subsys = {x:[] for x in set([x.subsystem.split("'")[1] for x in list(mat_model.reactions)])}
for r in list(mat_model.reactions):
    if r.id in reaction_ids:
        id_ind = np.where(ft.columns==r.id)[0][0]
        model_subsys[r.subsystem.split("'")[1]].append(id_ind) #.append(r.id)

def subsys_coverage(pddf):
    subsys = {pddf.index[kk]:
                  {k: pddf.iloc[kk, v].sum() for k, v in model_subsys.items()} for kk in range(pddf.shape[0])}
    subsys_df = pd.DataFrame().from_dict(subsys)
    return subsys_df

adipose_subsys = subsys_coverage(ft.loc[ft_tissues=="Adipose",:])
breast_subsys = subsys_coverage(ft.loc[ft_tissues=="Breast",:])
kidney_subsys = subsys_coverage(ft.loc[ft_tissues=="Kidney",:])
liver_subsys = subsys_coverage(ft.loc[ft_tissues=="Liver",:])

#deviation from mean
adipose_mean = np.mean(adipose_subsys,axis=1)
adipose_subsys_sub = adipose_subsys.subtract(adipose_mean,axis=0).abs()
adipose_subsys_dev = adipose_subsys_sub.div(adipose_mean,axis=0)*100 #5 means == 0 * 53 cols ->265

#at least a 10% deviation from mean subsystem coverage
(adipose_subsys_dev >= 0.1).sum()

# subsys_full = pd.concat([adipose_subsys,breast_subsys,kidney_subsys,liver_subsys],
#                         axis=1,keys=["Adipose","Breast","Kidney","Liver"])
# subsys_mean = np.mean(subsys_full,axis=1)
# subsys_full2 = subsys_full.subtract(subsys_mean,axis=0) #.abs()
# np.mean(subsys_full2.abs(),axis=1) > 0.1

### Functional comparison (tasks) ##
# tasks_subsys = np.array([t[-1][-1]["subsystem"] for t in h1_tasks.iterrows()
#                          if "subsystem" in t[-1][-1].keys()])

###
#set maximum growth rate
#test growth on Ham's media->growth task (GR1, 27)
# kidney = ft.iloc[np.where(ft_tissues=="Kidney")[0],:]

kidney_csm = read_sbml_model(MODEL_PATH) #mudar para "with model as csm" loop
for i,rct in ft.iloc[0,:].items():
    if rct==False: kidney_csm.reactions.get_by_id(i).knock_out()

kidney_csm.objective = "biomass_human"
#slim_optimize e mais rapido, da so objective value
kidney_csm.optimize().objective_value #produz biomassa




######################################################
# unfeasible = np.array([-1]*len(reactions)).reshape(len(reactions),)
# fva_df = np.zeros((len(reactions),ft.shape[0]*2))
# pfba_df = np.zeros((len(reactions),ft.shape[0]))
#
# for m in range(ft.shape[0]):
#     # model_b2 = read_sbml_model(MODEL_PATH)
#     with model as model_b2:
#         model_b2.objective = "biomass_human" #!!!
#         for r in range(ft.shape[1]):
#             if (ft.iloc[m, r] == False) and (ft.columns[r] not in essential_ids):
#                 model_b2.reactions.get_by_id(ft.columns[r]).bounds = (0, 0)
#
#         fba = model_b2.optimize()
#         if fba.status == "optimal": #and fba.objective_value > 1e-6
#             try:
#                 results_pfba = np.array(pfba(model_b2).fluxes)
#                 pfba_df[:, m] = results_pfba
#             except: pfba_df[:, m] = unfeasible
#
#             try: results_fva = np.array(flux_variability_analysis(model_b2, fraction_of_optimum=0.1))
#             except: results_fva = np.concatenate((unfeasible, unfeasible), axis=1)
#
#             fva_df[:, (m * 2):(m * 2) + 2] = results_fva
#         else:
#             pfba_df[:, m] = unfeasible
#             fva_df[:, (m * 2):(m * 2) + 2] = np.concatenate((unfeasible, unfeasible), axis=1)
#
# b20_FVA_path = os.path.join(ROOT_FOLDER,"DATA","results","b20_fva.csv")
# b20_pfba_path = os.path.join(ROOT_FOLDER,"DATA","results","b20_pfba.csv")
# pfba_df = pd.DataFrame(pfba_df)
# fva_df = pd.DataFrame(fva_df)
# pfba_df.to_csv(b20_pfba_path)
# fva_df.to_csv(b20_FVA_path)

# fva = pd.read_csv(FVA_path)
# max_zero, min_zero = abs(fva["1"]) < 1e-6, abs(fva["0"]) < 1e-6
# max_fwd, min_fwd = fva["1"] > 1e-6, fva["0"] > 1e-6
# max_rev, min_rev = fva["1"] < -1e-6, fva["0"] < -1e-6
# blk_df = max_zero & min_zero #0
# ess_df = (max_fwd & min_fwd) | (max_rev & min_rev) #88
# ((fva["0"].abs() > 1e-6) & (fva["1"].abs() > 1e-6 )).sum()

