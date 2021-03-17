# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 12:16:19 2020
@author: liamo
"""
import os, sys
from cobra.io import read_sbml_model, write_sbml_model
from cobra.flux_analysis import find_blocked_reactions
import pandas as pd
import numpy as np
import time

#windows
# ROOT_FOLDER = "C:\\Users\\liamo\\Documents\\BIOINF\\PRETHESIS"
# sys.path.append(os.path.join("C:\\Users\\liamo\\Documents\\BIOINF\\PRETHESIS\\troppo-dev\\src"))
# sys.path.append(os.path.join("C:\\Users\\liamo\\Documents\\BIOINF\\PRETHESIS\\cobamp-dev\\src"))
# HUMAN1_TASKS_PATH = os.path.join(ROOT_FOLDER,"Human-GEM","data","metabolicTasks", 'metabolicTasks_Essential.xlsx')
# HUMAN1_TASK_RESULTS_PATH = os.path.join(ROOT_FOLDER,"DATA","results", 'h1feb21og_ess_tasks.json')
# MODEL_PATH = os.path.join(ROOT_FOLDER,"Human-GEM","model", "Human-GEM_ct_feb21.xml")

ROOT_FOLDER = "/home/mpessoa/LINUX"
# MODEL_PATH = os.path.join(ROOT_FOLDER,"Human-GEM","model", "Human-GEM.xml")
MODEL_PATH = os.path.join(ROOT_FOLDER,"Human-GEM","model", "Human-GEM_ct.xml")
HUMAN1_TASKS_PATH = os.path.join(ROOT_FOLDER,"Human-GEM","data","metabolicTasks", 'metabolicTasks_Essential.xlsx')
HUMAN1_TASK_RESULTS_PATH = os.path.join(ROOT_FOLDER,"results", 'h1v1.5_tasks.json')
HUMAN1_TASK_RESULTS_CT_PATH = os.path.join(ROOT_FOLDER,"results", 'h1v1.5ct_tasks.json')

SOURCES_TO_ADD = "/home/mpessoa/LINUX/troppo-dev/src:/home/mpessoa/LINUX/troppo-dev:/home/mpessoa/LINUX/cobamp-dev/src:/home/mpessoa/LINUX/cobamp-dev"
for source in SOURCES_TO_ADD.split(':'):
    print('Adding source-code folder:', source)
    sys.path.append(source)

from cobamp.utilities.parallel import batch_run
from cobamp.wrappers.core import ConstraintBasedModelSimulator
from cobamp.wrappers.cobamp import *
from troppo.tasks.core import TaskEvaluator
from troppo.tasks.task_io import JSONTaskIO,ExcelTaskIO
from json import JSONEncoder

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

task_model = read_sbml_model(MODEL_PATH)
print("Reading model...")
task_list = get_human1_essential_tasks(task_model,HUMAN1_TASKS_PATH)
# create a task evaluator instance with the context specific model and the supplied task list and solver
task_eval = TaskEvaluator(model=task_model, tasks=task_list, solver='CPLEX')
# get task names (for future reference)
task_names = task_eval.tasks
# dir(task_eval.model),task_eval.current_task()

task_model.reactions.get_by_id('HMR_10024').bounds = (0, 1000)
task_model.objective = "biomass_human"
for task in task_list: task.mandatory_activity = []

print("Task evaluation...")
#template model task evaluation
#if a task does not pass in the original, it will not pass in the CSMs
with task_model:
    # tasks should be evaluated without open boundary reactions. we can easily close them on the COBRA model
    for k in task_model.boundary: k.knock_out()
    # use the batch_function from the TaskEvaluator class (takes the name of a loaded task, a params
    # dictionary with the task evaluator associated to the 'tev' key) and set the amount of threads to be used
    batch_res_tasks = batch_run(TaskEvaluator.batch_function, task_names, {'tev': task_eval}, threads=8)

# each element in the list of results in batch_res_tasks is a tuple of length 3 with the following:
# 0 - boolean flag representing the task evaluation
# 1 - Solution instance used to evaluate the task
# 2 - A dictionary with reactions supposed to be active mapped to True/False according to that criterion

# keep only items 0 and 2 of the task result - we don't need the flux distribution
task_eval_results = {k: (v[0], v[2]) for k, v in dict(zip(task_names, batch_res_tasks)).items()}

# save these results for later analysis as a JSON file
# with open(HUMAN1_TASK_RESULTS_PATH, 'w') as f:
#     f.write(JSONEncoder().encode(list(task_eval_results.items())))
with open(HUMAN1_TASK_RESULTS_CT_PATH, 'w') as f:
    f.write(JSONEncoder().encode(list(task_eval_results.items())))

print("Reactions required by each task...")
# reactions needed to pass the tasks
# similar to fva, obtains min and max flux for each reaction
# essential reactions are considered required
lentasks = len(task_eval.tasks)
lenrcts = len(task_eval.model.reaction_names)

cbms_fluxes_min_df = pd.DataFrame(np.zeros((lenrcts,lentasks)),index=task_eval.model.reaction_names, columns=task_eval.tasks)
cbms_status_min_df = pd.DataFrame(np.zeros((lenrcts,lentasks)),index=task_eval.model.reaction_names, columns=task_eval.tasks)

cbms_fluxes_max_df = pd.DataFrame(np.zeros((lenrcts,lentasks)),index=task_eval.model.reaction_names, columns=task_eval.tasks)
cbms_status_max_df = pd.DataFrame(np.zeros((lenrcts,lentasks)),index=task_eval.model.reaction_names, columns=task_eval.tasks)

# with task_model:
# tasks should be evaluated without open boundary reactions. we can easily close them on the COBRA model
for k in task_model.boundary: k.knock_out()

for i,t in enumerate(task_eval.tasks):
    cbms = ConstraintBasedModelSimulator(model=task_eval.model, simulation_function=cobamp_simulate,
                                         result_function=cobamp_simulation_result_function)
    task_eval.current_task = t #task_eval.tasks[i]
    print(f"Task {i+1} out of {lentasks} tasks")

    try:
        # c1_time = time.time()
        res_min = cbms.batch_simulate(func=cobamp_fba, bound_changes=[{}]*lenrcts,
                            objective_coefficients=[{k: 1} for k in task_eval.model.reaction_names], #task_eval.model.reaction_names
                            minimize=[True]*lenrcts, func_args=None, mp_threads=8)
        # c2_time = time.time()
        # print('Multi-threaded:', round(c2_time - c1_time,2), 'seconds',round((c2_time - c1_time)/60,2),"minutes")

        cbms_status_min_df[t] = [x[0] for x in res_min]
        cbms_fluxes_min_df[t] = [x[1] for x in res_min]

        res_max = cbms.batch_simulate(func=cobamp_fba, bound_changes=[{}]*lenrcts,
                            objective_coefficients=[{k: 1} for k in task_eval.model.reaction_names],
                            minimize=[False]*lenrcts, func_args=None, mp_threads=8)

        cbms_status_max_df[t] = [x[0] for x in res_max]
        cbms_fluxes_max_df[t] = [x[1] for x in res_max]

    except Exception as e:
      print(e)


cbms_out_fluxes_min = os.path.join(ROOT_FOLDER,"results","task_fva_test_fluxes_min.csv")
cbms_out_status_min = os.path.join(ROOT_FOLDER,"results","task_fva_test_status_min.csv")

cbms_out_fluxes_max = os.path.join(ROOT_FOLDER,"results","task_fva_test_fluxes_max.csv")
cbms_out_status_max = os.path.join(ROOT_FOLDER,"results","task_fva_test_status_max.csv")

cbms_fluxes_min_df.to_csv(cbms_out_fluxes_min)
cbms_status_min_df.to_csv(cbms_out_status_min)

cbms_fluxes_max_df.to_csv(cbms_out_fluxes_max)
cbms_status_max_df.to_csv(cbms_out_status_max)

# print((cbms_out_status_min==False).any().sum()) #0
# del cbms_out_status_min
# print((cbms_out_status_max==False).any().sum()) #2
# st_false = np.where(cbms_out_status_max==False)
# for rct in range(len(st_false[0])):
#     # print("Before:",stat_max.iloc[st_false[0][rct],st_false[1][rct]])
#     cbms_out_status_max.iloc[st_false[0][rct],st_false[1][rct]] = 0
#
# max_zero, min_zero = abs(cbms_out_fluxes_min) < 1e-6, abs(cbms_out_fluxes_max) < 1e-6
# max_fwd, min_fwd = cbms_out_fluxes_min > 1e-6, cbms_out_fluxes_max > 1e-6
# max_rev, min_rev = cbms_out_fluxes_min < -1e-6, cbms_out_fluxes_max < -1e-6
# blk_df = max_zero & min_zero
# ess_df = (min_fwd & max_fwd) | (min_rev & max_rev)
#
# # ess_df.loc[(ess_df==True).any(axis=1),:].index #170
# PROTECTED_EXTRA_PATH = os.path.join(ROOT_FOLDER,"task_rcts_protected.csv")
# protec_extra = pd.DataFrame(list(ess_df.loc[(ess_df==True).any(axis=1),:].index),columns=["ID"])
# protec_extra.to_csv(PROTECTED_EXTRA_PATH,index=False)
