import os, sys
import pandas as pd
import numpy as np
import re
import resource

if __name__ == '__main__':
    arg = [argval.split('=')[1] for argval in sys.argv[1:] if '-config=' in argval]
    # mp = '-no-mp' not in sys.argv[1:]
    mp = True

    if len(arg) > 0:
        ARGS_PATH = arg[0]
        with open(ARGS_PATH, 'r') as f:
            config_dict = {j[0].strip(): j[1].strip() for j in [k.split('=') for k in f.readlines()]}
    else:
        root_folder = "/home/mpessoa/LINUX"
        config_dict = {
            'ROOT_FOLDER': "/home/mpessoa/LINUX",
            'MODEL_PATH': os.path.join(root_folder,'Human-GEM_ct_nodrugex.xml'),
            # 'MODEL_PATH': os.path.join(root_folder,"Human-GEM","model", 'Human-GEM.xml'),
                    
            'DATA_PATH': os.path.join(root_folder,"DATA", 'gtexv8_more_rin6_gtl50_gtu90_lt50.csv'), #TPM
            # 'DATA_PATH': os.path.join(root_folder,"DATA", 'gtexv8_4tissues_rin6_gtl50_gtu90_lt50.csv'), 
            # 'DATA_PATH': os.path.join(root_folder,"DATA", 'gtexv8_median_by_TSAn_gtl50_gtu90_lt50.csv'), #GRP
            # 'DATA_PATH': os.path.join(root_folder,"DATA", 'gtexv8_median_gtl50_gtu90_lt50.csv'), #MED
            # 'DATA_PATH': os.path.join(root_folder, "DATA", 'gtexv8_4tissues_rin6_NOvarthresh_gtl50_gtu90_lt50.csv'), #TPM sem VAR thresh
            # 'DATA_PATH': os.path.join(root_folder,"DATA","gtexv8_4tissues_rin6_tpm.csv"), #test 1 TPM thres

            "METADATA_PATH": os.path.join(root_folder,"DATA", "metadata_clean.csv"),
            "NEW_METADATA_PATH": os.path.join(root_folder,"DATA", "allrin6_metadata_clean.csv"),
            
            'CS_MODEL_DF_FOLDER': os.path.join(root_folder,"results","CSMS"),
            'CS_MODEL_NAMES': 'gtl50_gtu90_lt50_ft_4tissue_novar_1tpmthresh_csm_nodrugex',
            'GENE_NAME_GRAB_PATTERN': 'ENSG[0-9]*',
            'PROTECTED_REACTION_LIST': 'biomass_human,HMR_10023,HMR_10024',
            'PROTECTED_REACTION_PATH': os.path.join(root_folder,"DATA", "task_rcts_protected_mat.csv"),
            "TASKS_PATH": os.path.join(root_folder,"Human-GEM","data","metabolicTasks", 'metabolicTasks_Essential.xlsx'),
            "TASK_RESULTS_PATH": os.path.join(root_folder,"results","TASKS","gtl50_gtu90_lt50_ft_4tissue_novar_1tpmthresh_task_nodrugex.json"),
            'INDEX_COLUMNS': '0,1,2,3', #'0,1'
            'NTHREADS': '8',
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

    from cobamp.wrappers.external_wrappers import get_model_reader
    from cobamp.utilities.parallel import batch_run, cpu_count
    from troppo.methods_wrappers import ReconstructionWrapper
    from troppo.omics.core import OmicsMeasurementSet, OmicsContainer

    from troppo.methods_wrappers import integration_strategy_map
    from troppo.omics.integration import MINSUM, MINMAX

    from itertools import product

    from troppo.tasks.core import TaskEvaluator
    from troppo.tasks.task_io import JSONTaskIO,ExcelTaskIO
    from json import JSONEncoder

    if not os.path.exists(config_dict['CS_MODEL_DF_FOLDER']): os.makedirs(config_dict['CS_MODEL_DF_FOLDER'])

    def get_human1_model():
        path, _ = urlretrieve('https://github.com/SysBioChalmers/Human-GEM/raw/master/model/Human-GEM.xml')
        model = read_sbml_model(path)
        # model.remove_metabolites([m for m in model.metabolites if m.compartment == 'x'])
        blocked = find_blocked_reactions(model)
        model.remove_reactions(blocked)
        write_sbml_model(model,'Human-GEM_ct.xml')
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

    print('Reading data...')
    
    tissue_meta_allrin = pd.read_csv(config_dict["NEW_METADATA_PATH"],index_col=0,header=0)
    tissue_meta = pd.read_csv(config_dict["METADATA_PATH"],index_col=0,header=0)
    new_tissues = ["Adipose - Subcutaneous","Breast - Mammary Tissue","Kidney - Cortex","Liver",
                "Whole Blood","Stomach","Lung","Brain - Cortex","Muscle - Skeletal",
                "Pancreas","Colon - Transverse"]
    
    #TPM
    # exp_data = pd.read_csv(config_dict['DATA_PATH'],
    #header=0,skiprows=[0,1],sep="\t",index_col="Name").drop(["Description"],axis=1).transpose().loc[new_tissues,:]
    exp_data = pd.read_csv(config_dict['DATA_PATH'],header=0,index_col=0)  #samples x genes
    # =list(map(int, config_dict['INDEX_COLUMNS'].split(','))))
    
    #median TPM
    # exp_data = exp_data.loc[["Adipose - Subcutaneous","Breast - Mammary Tissue","Kidney - Cortex","Liver"],:]
    # exp_data = exp_data.loc[new_tissues,:]
    
    #grouped
    # exp_data = pd.read_csv(config_dict['DATA_PATH'],header=0,index_col=[0,1,2])  #samples x genes
    # exp_data = exp_data.loc(axis=0)[new_tissues,:]
    # exp_data = exp_data.loc(axis=0)[["Adipose - Subcutaneous","Breast - Mammary Tissue","Kidney - Cortex","Liver"],:]
    
    #filter samples
    # tissue_ids = []
    # tissue_ids_max = []
    # for t in new_tissues: #tissue_meta["Tissue"].unique()
        # meta_indx = tissue_meta_allrin.iloc[np.where(tissue_meta_allrin["Tissue"] == t)[0], :].index #tissue_meta
        # sp_max = np.where(exp_data.index.str.contains("|".join(list(meta_indx)), case=False, regex=True))[0][22:24] #[20:40]  #[:53] 
        # tissue_ids_max.extend(sp_max)
    
        # sp = np.where(exp_data.index.str.contains("|".join(list(meta_indx)), case=False, regex=True))[0][:10]
        # tissue_ids.extend(sp)

    # exp_data = exp_data.iloc[tissue_ids,:]
    # exp_data = exp_data.iloc[list(set(tissue_ids_max).difference(set(tissue_ids))), :]
    # exp_data = exp_data.iloc[tissue_ids_max, :]


    # Create an omics measurement set object with the dataframe components
    omics_mset = OmicsMeasurementSet(exp_data.index, exp_data.columns, exp_data.values)

    # create a regex pattern to identify ensembl names in a string
    # Keep the ensembl gene ID only
    if 'GENE_NAME_GRAB_PATTERN' in config_dict.keys():
        print('\t','Grabbing gene names from pattern')
        ensembl_patt = re.compile(config_dict['GENE_NAME_GRAB_PATTERN'])
        omics_mset.column_names = [ensembl_patt.findall(k)[0] for k in omics_mset.column_names]

    data_dicts = {'_'.join(map(str, k)) if isinstance(k, tuple) else k: v for k, v in omics_mset.data.T.to_dict().items()}
    print("\tSamples:",len(data_dicts.keys()),"\n\tGenes per sample:",len(data_dicts[list(data_dicts.keys())[0]]))
    del exp_data, omics_mset

    print("Reading model...")
    if 'MODEL_PATH' in config_dict.keys() != '':
        model = read_sbml_model(config_dict['MODEL_PATH'])
        model.reactions.get_by_id('HMR_10024').bounds = (0, 1000)
        # blocked = find_blocked_reactions(model)
        # if len(blocked) > 0:
        # model.remove_reactions(blocked)
        # write_sbml_model(model, 'Human-GEM_ct.xml')
    else: get_human1_model()

    #protected reactions
    #minimal
    # protected = ['biomass_human', "HMR_10023", "HMR_10024"]
    # protected = list(map(lambda x: x.strip(),config_dict['PROTECTED_REACTION_LIST'].split(','))) \
    #     if ('PROTECTED_REACTION_LIST' in config_dict.keys()) else []

    #essencial task reactions
    protected = list(pd.read_csv(config_dict["PROTECTED_REACTION_PATH"], header=0)["ID"])

    print("Protected reactions:", len(protected))

    #### parameter setup
    #['fastcore',"tinit"]

    params = {'algorithms':['fastcore'],
              'strategies': {
                  'tinit': integration_strategy_map['adjusted_score'](protected),
                  'fastcore': integration_strategy_map['default_core'](0, protected)}, #omics=1 or gene scores=0
              'functions':{'minmax': MINMAX}}

    if 'OVERRIDE_COMBINATIONS' in config_dict:
        model_df_params = pd.read_csv(config_dict['OVERRIDE_COMBINATIONS'], index_col=0)
        models_to_reconstruct = [tuple(k) for k in model_df_params.values.tolist()]
    else: models_to_reconstruct = list(product(params['algorithms'], data_dicts.keys(), params['functions'].keys()))

    print('Reconstructing',len(models_to_reconstruct),'models...')

    for k in params['algorithms']:
        batch = [j for j in models_to_reconstruct if j[0] == k]
        print(len(batch),'models to reconstruct with',k)

    # model = read_sbml_model(config_dict['MODEL_PATH']) \
    #     if 'MODEL_PATH' in config_dict.keys() != '' \
    #     else get_human1_model()

    rw = ReconstructionWrapper(model, ttg_ratio=9999)
    params['rw'] = rw

    NTHREADS = int(cpu_count() if 'NTHREADS' not in config_dict.keys() else config_dict['NTHREADS'])
    safe_threads = {'tinit': max(1, NTHREADS // 16), 'fastcore': NTHREADS}

    print('Using',NTHREADS,'threads.')
    print('Available CPUs:',cpu_count())

    def reconstruct_model(options, params):
        # print('\tResource statistics before reconstruction with',options[0],':', resource.getrusage(resource.RUSAGE_SELF))
        alg, d, func = options
        data_dict, aofunc = d, params['functions'][func]

        oc_sample = OmicsContainer(omicstype='transcriptomics', condition='x',
                                   data=data_dict, nomenclature='custom')

        try:
            return params['rw'].run_from_omics(omics_data=oc_sample, algorithm=alg, and_or_funcs=aofunc,
                              integration_strategy=params['strategies'][alg], solver='CPLEX', raise_errors=False)

        except Exception as e:
            print(e)
            return {r: False for r in rw.model_reader.r_ids}


    # print('Resource statistics before reconstruction:',resource.getrusage(resource.RUSAGE_SELF))
    reconstructions = {}

    for k, v in safe_threads.items():
        batch = [j for j in models_to_reconstruct if j[0] == k]

        if len(batch) > 0:
            if mp:
                alg, dd, intf = zip(*batch)
                ddicts = [data_dicts[i] for i in dd]
                batch_args = list(zip(alg, ddicts, intf))
                print('\tResource statistics before reconstruction with',k,':', resource.getrusage(resource.RUSAGE_SELF))
                reconstructions.update(dict(zip(batch,batch_run(reconstruct_model, batch_args, params, threads=v))))
            else:
                for algo, ddi, intfx in batch:
                    data_dict = data_dicts[ddi]
                    print('Key = ', (algo, ddi, intfx))
                    reconstructions[(algo, ddi, intfx)] = reconstruct_model((algo, data_dict, intfx), params)

    print('Writing models')
    pd.DataFrame.from_dict(reconstructions, orient='index').\
        to_csv(os.path.join(config_dict['CS_MODEL_DF_FOLDER'], config_dict['CS_MODEL_NAMES'] + '.csv'))

    # print(pd.DataFrame.from_dict(reconstructions, orient='index').sum(axis=1))

    print("Task evaluation...")
    # Task evaluation #
    # parse tasks from a previously existing JSON
    # task_list = JSONTaskIO().read_task(config_dict["TASKS_PATH"])
    task_list = get_human1_essential_tasks(model)

    for task in task_list:
        # task.inflow_dict = {k: v if k not in task.outflow_dict.keys() else [-1000, 1000] for k, v in
        #                     task.inflow_dict.items()}
        # task.outflow_dict = {k: v for k, v in task.outflow_dict.items() if k not in task.inflow_dict.items()}
        task.mandatory_activity = []

    # read the original model to avoid compatibility issues with the tasks (e.g. missing metabolites from the block)
    task_model = read_sbml_model(config_dict["MODEL_PATH"])
    # tasks should be evaluated without open boundary reactions. we can easily close them on the COBRA model
    for k in task_model.boundary: k.knock_out()
    # get the names of all reactions in the model - this will be useful further on
    all_reactions = set([r.id for r in task_model.reactions])

    # since we have multiple models, we need to evaluate the tasks for each model
    # first, we create a structure to hold all of these results - a dictionary
    task_eval_results = {}

    # for each k (tuple with algorithm and sample information) and result (dict with reaction presences)...
    for k, result in reconstructions.items():
        # using with statements to change the COBRA model temporarily
        # this is done to knock-out reaction not appearing in the FASTCORE result
        with task_model as context_specific_model:
            protected = set([k for k, v in result.items() if v])  # get reactions included in the sample-specific model
            to_remove = all_reactions - protected  # get reactions except the protected ones
            for rid in to_remove:
                context_specific_model.reactions.get_by_id(rid).knock_out()  # knock-out reactions not in the model

            # create a task evaluator instance with the context specific model and the supplied task list and solver
            task_eval = TaskEvaluator(model=context_specific_model, tasks=task_list, solver='CPLEX')

            # get task names (for future reference)
            task_names = task_eval.tasks

            # use the batch_function from the TaskEvaluator class (takes the name of a loaded task, a params
            # dictionary with the task evaluator associated to the 'tev' key) and set the amount of threads to be used
            batch_res_tasks = batch_run(TaskEvaluator.batch_function, task_names, {'tev': task_eval}, threads=8)
        # each element in the list of results in batch_res_tasks is a tuple of length 3 with the following:
        # 0 - boolean flag representing the task evaluation
        # 1 - Solution instance used to evaluate the task
        # 2 - A dictionary with reactions supposed to be active mapped to True/False according to that criterion

        # keep only items 0 and 2 of the task result - we don't need the flux distribution
        task_csm_res = {k: (v[0], v[2]) for k, v in dict(zip(task_names, batch_res_tasks)).items()}
        print(k,"with", len(protected),"reactions and",
              len([v for k, v in task_csm_res.items() if v[0]]), 'tasks completed.')
        # assign this dictionary to it's sample on the master results dictionary
        task_eval_results[k] = task_csm_res

    # save these results for later analysis as a JSON file
    with open(config_dict["TASK_RESULTS_PATH"], 'w') as f:
        f.write(JSONEncoder().encode(list(task_eval_results.items())))
