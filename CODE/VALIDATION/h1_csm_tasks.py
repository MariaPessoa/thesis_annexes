import os, sys
import pandas as pd
import numpy as np

def get_human1_model():
    path, _ = urlretrieve('https://github.com/SysBioChalmers/Human-GEM/raw/master/model/Human-GEM.xml')
    model = read_sbml_model(path)
    # model.remove_metabolites([m for m in model.metabolites if m.compartment == 'x'])
    blocked = find_blocked_reactions(model)
    model.remove_reactions(blocked)
    write_sbml_model(model, 'Human-GEM_ct.xml')
    model.reactions.get_by_id('HMR_10024').bounds = (0, 1000)
    return model

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

if __name__ == '__main__':
    root_folder = "/home/mpessoa/LINUX"
    config_dict = {
        'ROOT_FOLDER': "/home/mpessoa/LINUX",
        # 'MODEL_PATH': os.path.join(root_folder, "Human-GEM", "model", 'Human-GEM_ct_nodrugex.xml'),
        'OG_MODEL_PATH': os.path.join(root_folder, "Human-GEM", "model", 'Human-GEM.xml'),
        'CS_MODELS': os.path.join(root_folder, "tINIT_GTEx_outputs.mat"),
        "TASKS_PATH": os.path.join(root_folder,"Human-GEM","data","metabolicTasks", 'metabolicTasks_Essential.xlsx'),
        "TASKS_FULL_PATH": os.path.join(root_folder, "Human-GEM", "data", "metabolicTasks", 'metabolicTasks_Full_task_sheet.xlsx'),
        "TASK_RESULTS_OUTPATH": os.path.join(root_folder,"results","TASKS", "tINIT_GTEx_7tissues_esstasks.json"),
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

    from cobra.io import read_sbml_model, write_sbml_model
    from cobra.util.solver import solvers
    from cobra.core.configuration import Configuration
    Configuration.solver = solvers['cplex']

    from cobra.flux_analysis import find_blocked_reactions
    from urllib.request import urlretrieve
    from cobamp.utilities.parallel import batch_run
    from troppo.tasks.core import TaskEvaluator
    from troppo.tasks.task_io import JSONTaskIO,ExcelTaskIO
    from json import JSONEncoder

    from scipy.io import loadmat
    from cobra.io.mat import from_mat_struct
    from cobra.flux_analysis import pfba,loopless_solution

    # print("Reading models from .mat file...")
    # #hardcoded
    # h1_csms = loadmat(config_dict["CS_MODELS"])["INIT_output"]
    # name_ids = pd.DataFrame([x[0][0] for x in h1_csms[0, 0][0]])
    # adipose = from_mat_struct(h1_csms[0, 0][1][0][0], model_id="adipose tissue")
    # breast = from_mat_struct(h1_csms[0, 0][1][np.where(name_ids == "breast")[0][0]][0], model_id="breast")
    # kidney = from_mat_struct(h1_csms[0, 0][1][np.where(name_ids == "kidney")[0][0]][0], model_id="kidney")
    # liver = from_mat_struct(h1_csms[0, 0][1][np.where(name_ids == "liver")[0][0]][0], model_id="liver")
    
    new_tissues = ["pancreas","stomach","lung","brain","muscle","blood","colon"]
    all_models = []
    for t in sorted(new_tissues):
        model = from_mat_struct(h1_csms[0, 0][1][np.where(name_ids == t)[0][0]][0], model_id=t)
        all_models.append(model)

    print("Reading model...")
    task_eval_results = {}
    task_model = read_sbml_model(config_dict["OG_MODEL_PATH"])
    all_reactions = set([r.id for r in task_model.reactions])
    
    task_list = get_human1_essential_tasks(task_model)
    task_list2 = get_human1_full_tasks(task_model)
    for k in task_model.boundary: k.knock_out()
    for task in task_list:
        task.mandatory_activity = []

    print("Task evaluation...")
    for m in all_models: #[adipose,breast,kidney,liver]
        with task_model as context_specific_model:
            to_remove = set([x.id for x in context_specific_model.reactions]).difference([x.id for x in m.reactions])
            for r in to_remove:
                if r in all_reactions: #extra
                    context_specific_model.reactions.get_by_id(r).knock_out()
    
            task_eval = TaskEvaluator(model=context_specific_model, tasks=task_list, solver='CPLEX')
            task_names = task_eval.tasks
            batch_res_tasks = batch_run(TaskEvaluator.batch_function, task_names, {'tev': task_eval}, threads=8)
    
        task_csm_res = {k: (v[0], v[2]) for k, v in dict(zip(task_names, batch_res_tasks)).items()}
        task_eval_results[m.id] = task_csm_res
    
    with open(config_dict["TASK_RESULTS_OUTPATH"], 'w') as f:
        f.write(JSONEncoder().encode(list(task_eval_results.items())))

    print("Full task evaluation...")
    for m in all_models: #[adipose,breast,kidney,liver]
        with task_model as context_specific_model:
            to_remove = set([x.id for x in context_specific_model.reactions]).difference([x.id for x in m.reactions])
            for r in to_remove:
                if r in all_reactions: #extra
                    context_specific_model.reactions.get_by_id(r).knock_out()
    
            task_eval = TaskEvaluator(model=context_specific_model, tasks=task_list2, solver='CPLEX')
            task_names = task_eval.tasks
            batch_res_tasks = batch_run(TaskEvaluator.batch_function, task_names, {'tev': task_eval}, threads=8)
    
        task_csm_res = {k: (v[0], v[2]) for k, v in dict(zip(task_names, batch_res_tasks)).items()}
        task_eval_results[m.id] = task_csm_res
    
    with open("tINIT_GTEx_7tissues_fulltasks.json", 'w') as f:
        f.write(JSONEncoder().encode(list(task_eval_results.items())))


