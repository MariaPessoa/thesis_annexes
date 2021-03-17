# -*- coding: utf-8 -*-
"""
Created on Jan 20
@author: liamo
"""
import pandas as pd
import numpy as np
import os, sys

from cobra.io import read_sbml_model
from cobra.flux_analysis import find_essential_reactions
from cobra.flux_analysis import pfba, flux_variability_analysis
from cobra.util.array import create_stoichiometric_matrix
sys.path.append(os.path.join("C:\\Users\\liamo\\Documents\\BIOINF\\PRETHESIS\\troppo-dev\\src"))
sys.path.append(os.path.join("C:\\Users\\liamo\\Documents\\BIOINF\\PRETHESIS\\cobamp-dev\\src"))
from troppo.tasks.task_io import ExcelTaskIO

ROOT_FOLDER = "C:\\Users\\liamo\\Documents\\BIOINF\\PRETHESIS\\"
MODEL_PATH = os.path.join(ROOT_FOLDER, "Human-GEM","model","Human-GEM_ct.xml")
HUMAN1_TASKS_PATH = os.path.join(ROOT_FOLDER,"Human-GEM","data","metabolicTasks", 'metabolicTasks_Essential.xlsx')
HUMAN1_TASK_RESULTS_PATH = os.path.join(ROOT_FOLDER,"DATA","results","h1v1.5_tasks.json")
HUMAN1_TASK_CT_RESULTS_PATH = os.path.join(ROOT_FOLDER,"DATA","results","h1v1.5ct_tasks.json")

cbms_out_fluxes_min = os.path.join(ROOT_FOLDER,"DATA","results","task_fva_test_fluxes_min.csv") #
cbms_out_status_min = os.path.join(ROOT_FOLDER,"DATA","results","task_fva_test_status_min.csv")
cbms_out_fluxes_max = os.path.join(ROOT_FOLDER,"DATA","results","task_fva_test_fluxes_max.csv")
cbms_out_status_max = os.path.join(ROOT_FOLDER,"DATA","results","task_fva_test_status_max.csv")

fluxes_min = pd.read_csv(cbms_out_fluxes_min,header=0,index_col=0)
stat_min = pd.read_csv(cbms_out_status_min,header=0,index_col=0)
fluxes_max = pd.read_csv(cbms_out_fluxes_max,header=0,index_col=0)
stat_max = pd.read_csv(cbms_out_status_max,header=0,index_col=0)

(stat_min==False).any().sum()
st_false = np.where(stat_min==False)
for rct in range(len(st_false[0])): fluxes_min.iloc[st_false[0][rct],st_false[1][rct]] = 0

(stat_max==False).any().sum()
st_false = np.where(stat_max==False)
for rct in range(len(st_false[0])):  fluxes_max.iloc[st_false[0][rct],st_false[1][rct]] = 0

max_zero, min_zero = abs(fluxes_min) < 1e-6, abs(fluxes_max) < 1e-6
max_fwd, min_fwd = fluxes_min > 1e-6, fluxes_max > 1e-6
max_rev, min_rev = fluxes_min < -1e-6, fluxes_max < -1e-6
blk_df = max_zero & min_zero
ess_df = (min_fwd & max_fwd) | (min_rev & max_rev)

ess_df.loc[(ess_df==True).any(axis=1),:].index #feb21 ->170, 1.5 -> 87
PROTECTED_EXTRA_PATH = os.path.join(ROOT_FOLDER,"DATA","task_rcts_protected.csv")
protec_extra = pd.DataFrame(list(ess_df.loc[(ess_df==True).any(axis=1),:].index),columns=["ID"])
protec_extra.to_csv(PROTECTED_EXTRA_PATH)


model = read_sbml_model(MODEL_PATH)
model.objective = "biomass_human"
model.reactions.get_by_id('HMR_10024').bounds = (0, 1000)
model.reactions.get_by_id('biomass_human').bounds = (0, 1000)

def get_human1_essential_tasks(model,taskfile):
    # URL = 'https://github.com/SysBioChalmers/Human-GEM/raw/master/data/metabolicTasks/metabolicTasks_Essential.xlsx'
    # path, _ = urlretrieve(URL)
    task_list = ExcelTaskIO().read_task(taskfile)

    metab_map = {m.name + '[' + m.compartment + ']': m.id if m.compartment != 'x' else m.id[:-1] + 's' for m in
                 model.metabolites}
    metab_map['NEFA blood pool in[x]'] = 'm02560s'
    replace_func = lambda x: metab_map[x] if x in metab_map.keys() else x

    for t in task_list: t.id_replace(replace_func)  # only works for inflow_dict/outflow_dict
    for t in task_list:
        t.reaction_dict = {k: [{metab_map[i]: j for i, j in v[0].items() if i in metab_map.keys()}, v[1]]
                           for k, v in t.reaction_dict.items()}

    task_list[27].outflow_dict.update(task_list[27].inflow_dict)
    del task_list[27].outflow_dict['ALLMETSIN[s]']
    return task_list

h1_tasks = get_human1_essential_tasks(model,HUMAN1_TASKS_PATH)
#h1feb21_ct falha BS15 / h1_tasks[6]
h1_task_res = pd.read_json(HUMAN1_TASK_RESULTS_PATH)
h1_task_ct_res = pd.read_json(HUMAN1_TASK_CT_RESULTS_PATH)

tasks_met_out = []
for row in range(len(h1_tasks)):
    tasks_met_out.extend(list(h1_tasks[row].outflow_dict.keys()))
tasks_met_out2 = set(tasks_met_out) #len = 97, linux doesnt like it the normal way??

modelS = create_stoichiometric_matrix(model,"DataFrame") #metabolites x reactions
#only reactions involving the metabolites of interest
modelS_tasks = modelS.loc[tasks_met_out2,:]
modelS_tasks.loc[(modelS_tasks>=1).any(axis=1),:].shape #2 dont have a production reaction
modelS_tasks.loc[~(modelS_tasks>=1).any(axis=1),:].index

#task evaluation blocks boundary reactions
bdry = [r.id for r in model.boundary]
modelS_tasks = modelS_tasks.loc[:,~modelS.columns.isin(set(bdry))]
#reactions that produce a target metabolite?
modelS_tasks = modelS_tasks.loc[:,(modelS_tasks>=1).any()]


