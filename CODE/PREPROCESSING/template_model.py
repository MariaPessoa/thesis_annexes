# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 12:16:19 2020
@author: liamo
"""
import os, sys
import pandas as pd
import numpy as np

root_folder = "/home/mpessoa/LINUX"
config_dict = {
    'root_folder': "/home/mpessoa/LINUX",
    # 'MODEL_PATH': os.path.join(root_folder, 'Human-GEM_ct_nodrugex.xml'),
    'MODEL_PATH': os.path.join(root_folder, "DATA","UHLEN", 'adipose.xml'),
    'OG_MODEL_PATH': os.path.join(root_folder, 'Human-GEM_ct.xml'),
    "TASKS_PATH": os.path.join(root_folder, "Human-GEM", "data", "metabolicTasks", 'metabolicTasks_Essential.xlsx'),
    "TASKS_FULL_PATH": os.path.join(root_folder, "Human-GEM", "data", "metabolicTasks", 'metabolicTasks_Full_task_sheet.xlsx'),
    # "TASK_RESULTS_OUTPATH": os.path.join(root_folder, "results", "TASKS", "Uhlen_adipose_full.json"),
    "TASK_RESULTS_OUTPATH": os.path.join(root_folder, "DATA","UHLEN", "Uhlen_adipose_full.json"),
    'SOURCES_TO_ADD': "/home/mpessoa/LINUX/troppo-dev/src:/home/mpessoa/LINUX/troppo-dev:/home/mpessoa/LINUX/cobamp-dev/src:/home/mpessoa/LINUX/cobamp-dev"
}

print('Configuration:')
for k, v in config_dict.items():
    print("\t" + k, '=', v)

if 'SOURCES_TO_ADD' in config_dict:
    sources_to_add = config_dict['SOURCES_TO_ADD'].split(':')
    for source in sources_to_add:
        print('Adding source-code folder:', source)
    sys.path.extend(sources_to_add)


from itertools import chain
from cobra.flux_analysis import find_blocked_reactions
from cobra.io import read_sbml_model, write_sbml_model, load_matlab_model
from cobra.util.solver import solvers
from cobra.core.configuration import Configuration
Configuration.solver = solvers['cplex']

from cobamp.utilities.parallel import batch_run
from troppo.tasks.core import TaskEvaluator
from troppo.tasks.task_io import JSONTaskIO, ExcelTaskIO
from json import JSONEncoder
from troppo.tasks.core import Task

print("Reading base model...")
model = read_sbml_model(config_dict["OG_MODEL_PATH"])
all_reactions = [x.id for x in list(model.reactions)]

#gerar o ct_nodrugex
print("Saving template model...")
model_xls = pd.ExcelFile(os.path.join(root_folder,"Human-GEM","model", "Human-GEM.xlsx")).parse(sheet_name='RXNS').drop('#', axis=1)
subsystems = {s:model_xls.loc[model_xls['SUBSYSTEM']==s]['ID'].tolist() for s in model_xls['SUBSYSTEM'].unique().tolist()}

meta_reacs = {sub:set(chain(
    *[set(m.id for m in model.reactions.get_by_id(r).metabolites.keys()) for r in reacs if r in all_reactions])
    ) for sub,reacs in subsystems.items()}
drugs_from_exchange = meta_reacs['Drug metabolism'] & meta_reacs['Exchange/demand reactions'] - \
                      set(chain(*[v for k,v in meta_reacs.items() if k not in ['Drug metabolism', 'Exchange/demand reactions']]))
reactions_drugs_from_exchange = [r.id for r in model.boundary if len(set([m.id for m in r.metabolites]) & drugs_from_exchange) > 0]
# len(reactions_drugs_from_exchange) #146

for r in reactions_drugs_from_exchange:
    model.reactions.get_by_id(r).bounds = (0,0)
model.objective = "biomass_human"
blocked = find_blocked_reactions(model)
print("blocked:",len(blocked)) #516
model.remove_reactions(blocked)
write_sbml_model(model,'Human-GEM_ct_nodrugex.xml')

def get_human1_essential_tasks(model):
    # URL = 'https://github.com/SysBioChalmers/Human-GEM/raw/master/data/metabolicTasks/metabolicTasks_Essential.xlsx'
    # path, _ = urlretrieve(URL)
    task_list = ExcelTaskIO().read_task(config_dict["TASKS_PATH"])

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


def get_human1_full_tasks(model):
    # URL = 'https://github.com/SysBioChalmers/Human-GEM/raw/master/data/metabolicTasks/metabolicTasks_Essential.xlsx'
    # path, _ = urlretrieve(URL)
    task_list = ExcelTaskIO().read_task(config_dict["TASKS_FULL_PATH"])

    metab_map = {m.name + '[' + m.compartment + ']': m.id if m.compartment != 'x' else m.id[:-1] + 's' for m in
                 model.metabolites}
    replace_func = lambda x: metab_map[x] if x in metab_map.keys() else x

    for t in task_list: t.id_replace(replace_func)  # only works for inflow_dict/outflow_dict
    for t in task_list:
        t.reaction_dict = {k: [{metab_map[i]: j for i, j in v[0].items() if i in metab_map.keys()}, v[1]]
                           for k, v in t.reaction_dict.items()}
    return task_list
    

print("Reading model...")
task_model = read_sbml_model(config_dict["MODEL_PATH"]) #config_dict["OG_MODEL_PATH"]

print("Task Evaluation...")
# task_list = get_human1_essential_tasks(model)
# task_list = get_human1_full_tasks(task_model)

task_list = ExcelTaskIO().read_task(config_dict["TASKS_FULL_PATH"])

for task in task_list:
    # task.inflow_dict = {k: v if k not in task.outflow_dict.keys() else [-1000, 1000] for k, v in
    #                     task.inflow_dict.items()}
    # task.outflow_dict = {k: v for k, v in task.outflow_dict.items() if k not in task.inflow_dict.items()}
    task.mandatory_activity = []

# read the original model to avoid compatibility issues with the tasks (e.g. missing metabolites from the block)
# tasks should be evaluated without open boundary reactions. we can easily close them on the COBRA model
for k in task_model.boundary: k.knock_out()
all_reactions = set([r.id for r in task_model.reactions])
task_eval = TaskEvaluator(model=task_model, tasks=task_list, solver='CPLEX')
task_names = task_eval.tasks

# use the batch_function from the TaskEvaluator class (takes the name of a loaded task, a params
# dictionary with the task evaluator associated to the 'tev' key) and set the amount of threads to be used
# try:
batch_res_tasks = batch_run(TaskEvaluator.batch_function, task_names, {'tev': task_eval}, threads=8)
# except Exception as e: #se der erro depois o zip nao da pq lens diff
#     print(e)
# keep only items 0 and 2 of the task result - we don't need the flux distribution
task_eval_results = {k: (v[0], v[2]) for k, v in dict(zip(task_names, batch_res_tasks)).items()}

with open(config_dict["TASK_RESULTS_OUTPATH"], 'w') as f:
    f.write(JSONEncoder().encode(list(task_eval_results.items())))



