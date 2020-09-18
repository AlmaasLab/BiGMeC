#!/usr/bin/env python
# coding: utf-8
"""
Authors:
  - Snorre Sulheim, snorres.sulheim@sintef.no

Date: 17.09.2020
Lisence: CC-BY-4.0

This file contains functions used to incorporate metabolic pathways and run optknock to predict optimal strain design
"""

import cobra
import cameo
from pathlib import Path
import pandas as pd
from cameo.strain_design import OptKnock
import time
import numpy as np
import logging
import datetime


SOLVER = "gurobi"
OPTKNOCK_BIOMASS_FRACTION = 0.999
MINIMUM_FLUX = 0.01
EPSILON = 1e-2

def add_metabolic_pathway(model, pathway_fn):
    pathway = cobra.io.load_json_model(pathway_fn)
    model.merge(pathway)
    return pathway.description


def BGC_optknock(model, pathway_fn, optknock_reactions = [], fraction_of_optimum = 0.9, biomass_rxn_id = "BIOMASS_SCO_tRNA", target_id = "DM_secondary_metabolite", brute_force = True):
    model.solver = SOLVER
    bgc_type = add_metabolic_pathway(model, pathway_fn)
    exclude_reactions = [r.id for r in model.reactions if not r.id in optknock_reactions]

    model.objective = model.reactions.get_by_id(biomass_rxn_id)
    max_growth = model.slim_optimize()
    model.reactions.get_by_id(biomass_rxn_id).lower_bound = max_growth*fraction_of_optimum

    if brute_force:
        results = single_optknock(model, optknock_reactions, biomass_rxn_id, target_id, max_growth*fraction_of_optimum)
        results = pd.DataFrame(results)
        results.columns = ["ID", "Growth", "Production"]
        results.set_index("ID", inplace = True)
    
        improved_results = results.loc[results["Production"] > results.loc["WT", "Production"] + EPSILON, :]
        improved_results = improved_results/results.loc["WT", :]
        improved_results["BGC type"] = bgc_type
        #result = results.sort_values("Production", ascending=False).iloc[:10,:]/results.
        #result.append(results.loc[results["ID"]=="WT", :])

        return improved_results.sort_values("Production", ascending = False)
    else:
        try:
            optknock = OptKnock(model, exclude_reactions = exclude_reactions, fraction_of_optimum=None)
            result = optknock.run(max_knockouts=1, target=target_id, biomass=biomass_rxn_id, max_results = 1)
        except KeyError:
            logging.warning("No secondary metabolite in model", pathway_fn)
        return result


def folder_BGC_optknock(model_fn, folder, fraction_of_optimum = 0.9, 
                        biomass_rxn_id = "BIOMASS_SCO_tRNA", target_id = "DM_secondary_metabolite", 
                        results_folder = "."):
    model = cobra.io.read_sbml_model(model_fn)
    max_growth = model.slim_optimize()
    optknock_reactions = optknock_prepare(model, max_growth*fraction_of_optimum)
    #optknock_reactions = optknock_reactions[::10]

    folder = Path(folder)
    all_results = []
    for i, pathway_fn in enumerate(folder.glob("*.json")):
        logging.info(pathway_fn)
        with model:
            result = BGC_optknock(model, str(pathway_fn), optknock_reactions = optknock_reactions, 
                                  fraction_of_optimum = fraction_of_optimum, brute_force = True)
            save_fn = results_folder + "/brute_optknock_{0}_{1}.csv".format(str(pathway_fn.stem), fraction_of_optimum)
            result.to_csv(save_fn)
            logging.info(pathway_fn, result)
            result["BGC"] = str(pathway_fn.stem)
            all_results.append(result)

    df = pd.concat(all_results)
    logging.info(df)
    df.to_csv("{0}/all_optknock_results_{1}.csv".format(results_folder, fraction_of_optimum))

def single_optknock(model, reaction_list, biomass_rxn_id, target_reaction, minimum_growth = 0.01):
    single_optknock_list = []

    # Wild type
    growth, production = optknock_brute_force(model, biomass_rxn_id, target_reaction, minimum_growth)
    single_optknock_list.append(["WT", growth, production])
    print("Wild-type: ", growth, production)

    for r_id in reaction_list:
        with model as model:
            model.reactions.get_by_id(r_id).knock_out()
            growth, production = optknock_brute_force(model, biomass_rxn_id, target_reaction, minimum_growth)
        single_optknock_list.append([r_id, growth, production])
    return single_optknock_list

def optknock_brute_force(model, biomass_rxn_id, target_reaction, minimum_growth):
    model.objective = model.reactions.get_by_id(biomass_rxn_id)
    growth = model.slim_optimize()
    #print("Growth: ", growth)
    if growth < minimum_growth:
        print("Zero growth")
        return 0, 0
    elif np.isnan(growth):
        print("No solution for model")
        return 0, 0
    else:
        with model:
            model.reactions.get_by_id(biomass_rxn_id).lower_bound = growth * OPTKNOCK_BIOMASS_FRACTION
            model.objective = model.reactions.get_by_id(target_reaction)
            production = model.slim_optimize()
        return growth, production


def optknock_prepare(model, minimum_growth = MINIMUM_FLUX):
    # Only keep reactions with gene associations
    gene_reaction_ids = []
    for reaction in model.reactions:
        if len(reaction.gene_reaction_rule) and reaction.gene_reaction_rule != "s0001":
            gene_reaction_ids.append(reaction.id)

    # Remove blocked reactions
    blocked_reactions = cobra.flux_analysis.find_blocked_reactions(model)
    optknock_reactions = [r for r in gene_reaction_ids if not r in blocked_reactions]
        
    # Get list of essential_reactions
    essential_reactions = get_essential_reactions(model, reaction_list = optknock_reactions, minimum_growth = minimum_growth)

    #Remove essential reactions
    optknock_reactions = [r for r in optknock_reactions if not r in essential_reactions]

    #Remove exchanges
    exchanges = [r.id for r in model.exchanges]
    optknock_reactions = [r for r in optknock_reactions if not r in exchanges]
    return optknock_reactions


# def optknock_brute_force(model, target_reaction, solver  = "gurobi", n_jobs = -1, output = "optknock_brute_force.csv", minimum_growth = 0.01):
    
#     s0 = time.time()
#     s_arr = []

#         # model.solver = solver
#     if n_jobs == -1:
#         n_jobs = multiprocessing.cpu_count()

#     print("Number of jobs; ", n_jobs)
    

#     # optknock_reactions = optknock_reactions[::20]
#     print("Number of optknock_reactions:", len(optknock_reactions))
#     # optknock_triplets = list(itertools.combinations(np.arange(len(optknock_reactions)), 3))

#     print(s_arr)

#     # Find lethal pairs
    
#     # print("Found lethal pairs: ", len(lethal_pairs), len(nonlethal_pairs))

#     triple_optknock_df = triple_optknock(model, optknock_reactions, target_reaction, n_jobs, minimum_growth)
#     s_arr.append(time.time()-s0)
#     print("Finished triplets")
    
#     lethal_pairs, nonlethal_pairs = get_lethal_and_nonlethal_pairs(model, optknock_reactions)
#     double_optknock_list = double_optknock(model, nonlethal_pairs, target_reaction, minimum_growth)
#     print("Finished doubles")
    
#     s_arr.append(time.time()-s0)
#     single_optknock_list = single_optknock(model, optknock_reactions, target_reaction, minimum_growth)
#     s_arr.append(time.time()-s0)
#     print(s_arr)

#     optknock_list = single_optknock_list + double_optknock_list
    

#     #Store data
#     df12 = pd.DataFrame(optknock_list, columns = ["Reaction 1", "Reaction 2", "Reaction 3", "Growth", "Production"])
#     # df.columns = ["Reaction 1", "Reaction 2", "Reaction 3", "Growth", "Production"]
#     df_growth = df12[df12["Growth"]>0]

#     df = pd.concat([triple_optknock_df, df_growth], ignore_index = True)


#     fn = "env_"+str(env)+"_"+output
#     path = os.path.join(os.getcwd(), fn)
#     df.to_csv(path, index = False)
#     print("Saved!")


def get_essential_reactions(model, reaction_list = None, minimum_growth = MINIMUM_FLUX):
    if isinstance(reaction_list, list):
        reaction_list = [model.reactions.get_by_id(x) for x in reaction_list]

    single_reaction_ko = cobra.flux_analysis.single_reaction_deletion(model, reaction_list = reaction_list)
    single_reaction_ko[~single_reaction_ko["growth"].isna()] = 0
    essential_reactions = list(single_reaction_ko[single_reaction_ko["growth"] < minimum_growth].index)
    return essential_reactions



if __name__ == '__main__':
    model_fn = "../Models/Sco-GEM.xml"
    pathway_fn = "../Data/constructed_pathways/1.json"
    folder = "../Data/constructed_pathways/"
    results_folder = "../Data/knockouts"


    time = datetime.datetime.now()

    logging.basicConfig(filename='strain_design_{0}_{1}.log'.format(time.month, time.day), filemode='w', format='%(name)s - %(levelname)s - %(message)s', level = logging.INFO)
    
    folder_BGC_optknock(model_fn, folder, fraction_of_optimum = 0.5, results_folder = results_folder)

