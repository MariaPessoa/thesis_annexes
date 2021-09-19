import os, sys
import pandas as pd
import numpy as np

if __name__ == '__main__':
    root_folder = "/home/mpessoa/LINUX"
    config_dict = {
        'ROOT_FOLDER': "/home/mpessoa/LINUX",
        'MODEL_PATH': os.path.join(root_folder, 'Human-GEM_ct_nodrugex.xml'),
        # 'OG_MODEL_PATH': os.path.join(root_folder, 'Human-GEM_ct.xml'),
        'OG_MODEL_PATH': os.path.join(root_folder,"Human-GEM","model", 'Human-GEM.xml'),
        "METADATA_PATH": os.path.join(root_folder,"DATA", "metadata_clean.csv"),
        'CS_MODEL_DF_FOLDER': os.path.join(root_folder,"results","CSMS"),
        'CS_MODEL_NAMES': ['gtl50_gtu90_lt50_ft_4tissues_rin6_prmat_csm_nodrugex','gtl50_gtu90_lt50_ft_extra_rin6_prmat_csm_nodrugex'],
        "TASKS_PATH": os.path.join(root_folder,"Human-GEM","data","metabolicTasks", 'metabolicTasks_Essential.xlsx'),
        "TASKS_FULL_PATH": os.path.join(root_folder, "Human-GEM", "data", "metabolicTasks", 'metabolicTasks_Full_task_sheet.xlsx'),
        "TASK_RESULTS_OUTPATH": os.path.join(root_folder,"results","TASKS","gtl50_gtu90_lt50_fttinit_med_1tpmthresh_full_nodrugex.json"),
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
    from troppo.tasks.core import Task


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


    ##
    # if not os.path.exists(config_dict['CS_MODEL_DF_FOLDER']): os.makedirs(config_dict['CS_MODEL_DF_FOLDER'])
    print('Reading data...')
    all_tissues = ['adipose tissue', 'brain','breast', 'blood', 'colon', 'kidney', 'liver', 'lung', 
       'muscle', 'pancreas', 'stomach']
    #"Unnamed: 0",
    csms = pd.read_csv(os.path.join(root_folder,"results","CSMS",
                          "gtl50_gtu90_lt50_fttinit_med_1tpmthresh_csm_nodrugex.csv"),
                          header=0,index_col=[0,1]).drop(["Unnamed: 2"], axis=1)

    # if len(config_dict["CS_MODEL_NAMES"]) == 1:
        # csms = pd.read_csv(os.path.join(config_dict["CS_MODEL_DF_FOLDER"], config_dict["CS_MODEL_NAMES"] + ".csv"),
                          # header=0, index_col=1).drop(["Unnamed: 0", "Unnamed: 2"], axis=1)
    # else: #can be empty, doesnt check
        # models = []
        # for m in config_dict["CS_MODEL_NAMES"]:
            # models.append(pd.read_csv(os.path.join(config_dict["CS_MODEL_DF_FOLDER"], m + ".csv"),
                          # header=0, index_col=1).drop(["Unnamed: 0", "Unnamed: 2"], axis=1))
            # #keys=list(range(len(config_dict["CS_MODEL_NAMES"])))?
            # csms = pd.concat(models,axis=0) #,keys=config_dict["CS_MODEL_NAMES"]
            # # csms = csms.reindex(csms_ft.index) #ensure same order as ft models
    
    print("csms df shape:",csms.shape)
    print("Reading model...")
    if 'MODEL_PATH' in config_dict.keys() != '': #"MODEL_PATH"
        task_model = read_sbml_model(config_dict['MODEL_PATH']) #MODEL_PATH
        task_model.reactions.get_by_id('HMR_10024').bounds = (0, 1000)
        # blocked = find_blocked_reactions(model)
        # model.remove_reactions(blocked)
        # write_sbml_model(model, 'Human-GEM_ct.xml')
    else: get_human1_model()
    print("Model reaction count (ignores blocked):",len(list(task_model.reactions)))

    print("Task evaluation...")
    # task_list = get_human1_essential_tasks(task_model)
    task_list = get_human1_full_tasks(task_model)

    for task in task_list:
        # task.inflow_dict = {k: v if k not in task.outflow_dict.keys() else [-1000, 1000] for k, v in
        #                     task.inflow_dict.items()}
        # task.outflow_dict = {k: v for k, v in task.outflow_dict.items() if k not in task.inflow_dict.items()}
        task.mandatory_activity = []

    # read the original model to avoid compatibility issues with the tasks (e.g. missing metabolites from the block)
    # tasks should be evaluated without open boundary reactions. we can easily close them on the COBRA model
    for k in task_model.boundary: k.knock_out()
    all_reactions = set([r.id for r in task_model.reactions])
    task_eval_results = {}

    #csms / csms_ft
    for i,row in csms.iterrows(): #expects reaction content
        # using with statements to change the COBRA model temporarily
        # this is done to knock-out reaction not appearing in the FASTCORE result
        with task_model as context_specific_model:
            to_remove = row[(row==False)].index.to_list()
            to_remove+=list(set(all_reactions).difference(csms.columns.to_list())) #H1_csms only
            for r in to_remove: 
                # if r in all_reactions:
                context_specific_model.reactions.get_by_id(r).knock_out() # knock-out reactions not in the model

            # create a task evaluator instance with the context specific model and the supplied task list and solver
            task_eval = TaskEvaluator(model=context_specific_model, tasks=task_list, solver='CPLEX')
            task_names = task_eval.tasks

            # use the batch_function from the TaskEvaluator class (takes the name of a loaded task, a params
            # dictionary with the task evaluator associated to the 'tev' key) and set the amount of threads to be used
            batch_res_tasks = batch_run(TaskEvaluator.batch_function, task_names, {'tev': task_eval}, threads=8)

        # keep only items 0 and 2 of the task result - we don't need the flux distribution
        task_csm_res = {k: (v[0], v[2]) for k, v in dict(zip(task_names, batch_res_tasks)).items()}
        task_eval_results[i] = task_csm_res

    with open(config_dict["TASK_RESULTS_OUTPATH"], 'w') as f:
        f.write(JSONEncoder().encode(list(task_eval_results.items())))

