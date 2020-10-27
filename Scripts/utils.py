#!/usr/bin/env python
# coding: utf-8
"""
Authors:
  - Snorre Sulheim, snorres.sulheim@sintef.no

Date: 17.09.2020
Lisence: CC-BY-4.0

Useful utility functions

"""

import pandas as pd
import cobra
from pathlib import Path
import json
import re

name_fix_dict = {
            "nadph": "NADPH",
            "nadh": "NADH",
            "nad": "NAD+",
            "nadp": "NADP+",
            "co2": "CO2",
            "coa": "CoA",            
            "ppi": "Diphosphate",
            "h2o": "H2O",
            "amp": "AMP",
            "atp": "ATP",
            "Generic amino acid": "AA-X",
            "mmcoa__R": "Methylmalonyl-CoA"
}
name_fix_dict_reverse = {value:key for (key, value) in name_fix_dict.items()}

def pathway_to_table(json_fn, save_folder = None, include_backbone = False):
    model = cobra.io.load_json_model(str(json_fn))
    backbone_list = []
    for r in model.reactions:
        split_string = r.id.split("_")
        if len(split_string) == 2:
            bgc_type, j = split_string
            try:
                int(j)
            except ValueError:
                # If j can't be converted to an integer this is not one of the backbone reactions
                continue

            sub_string, product_string = r.reaction.split("-->")
            # Remove e.g. NRPS_1 + .. from the reaction string
            if not include_backbone:
                search_string = r"{0}+\_\d+".format(bgc_type)
                sub_string = re.sub(search_string, "", sub_string).replace("+  +", "+").strip(" + ").strip()
                product_string   = re.sub(search_string, "", product_string).replace("+  +", "+").strip(" + ").strip()
            
            # Name

            product_mets = []
            substrate_mets = []
            for m, i in r.metabolites.items():
                if m.name != bgc_type:
                    if m.id[:-2] in name_fix_dict:
                        name = name_fix_dict[m.id[:-2]]
                    else:
                        name = m.name
                    if i < 0:
                        substrate_mets.append(name)
                    else:
                        product_mets.append(name)

            full_ID_string = "{0} -> {1}".format(sub_string, product_string)

            product_name_string = " + ".join(product_mets)
            substrate_names_string = " + ".join(substrate_mets)
            full_name_string = "{0} -> {1}".format(substrate_names_string, product_name_string)
            gene = r.gene_reaction_rule
                
            backbone_list.append([int(j), sub_string, product_string, full_ID_string,
                                  substrate_names_string, product_name_string, full_name_string, gene])
    df = pd.DataFrame(backbone_list, columns = ["N", "Substrate IDs", "Product IDs", "Full reaction (IDs)",
                    "Substrate names", "Product names","Full reaction (names)", "gene"])
    df.sort_values(by = "N", inplace  =True)
    df.set_index("N", inplace = True)
    if save_folder:
        name = str(Path(json_fn.stem))+".csv"
        fn = Path(save_folder) / name
        df.to_csv(fn)
    else:
        print(json_fn)
        print(df)

def all_pathways_to_table(json_folder, csv_folder):
    csv_folder = Path(csv_folder)
    csv_folder.mkdir(exist_ok = True)

    for json_fn in Path(json_folder).glob("*.json"):
        pathway_to_table(json_fn, csv_folder)

def fix_BiGG_json():
    """
    The universal bigg models downloaded from http://bigg.ucsd.edu/data_access has annotations as a list, but it has to be as dict.
    """
    model_fn = "../Models/BiGG_universal_model.json"
    
    # Read in json file
    with open(model_fn, "r") as f:
        data = json.load(f)

    # change annotations from list to dict
    for key,key_data in data.items():
        if isinstance(key_data, list):
            for element_dict in key_data:
                try:
                    annotations = element_dict["annotation"]
                except KeyError:
                    continue
                except TypeError:
                    continue
                else:
                    element_dict["annotation"] = {k:v for (k,v) in annotations}

    # Store json
    with open(model_fn, "w") as f:
        json.dump(data, f)

def make_pathway_from_csv():
    model_fn = '../Models/Sco-GEM.xml'
    ref_model = cobra.io.read_sbml_model(model_fn)
    
    fn = "../Data/pathway_comparison.xlsx"

    df_dict = pd.read_excel(fn, sheet_name = None, skiprows = 14)
    for key, df in df_dict.items():
        print("\n", "# ", key)

        m_sum_dict = {}
        for string in df["Real reaction"]:
            if isinstance(string, str) and len(string.lstrip().rstrip()):
                # print(string)
                m_dict = reaction_string_to_mets(string, ref_model)
                #print(m_dict)
                for m, value in m_dict.items():
                    try:
                        m_id = m.id
                    except AttributeError:
                        m_id = m
                    try:
                        current_value = m_sum_dict[m_id]
                    except KeyError:
                        m_sum_dict[m_id] = value
                    else:
                        m_sum_dict[m_id] = current_value + value
        for key, value in m_sum_dict.items():
            print("{0}: {1}".format(key, value))


def reaction_string_to_mets(string, model):
    substrate_string, product_string = string.split("->")
    substrates = substrate_string.split(" + ")
    products = product_string.split(" + ")

    m_dict = {}
    for string, k in zip([substrates, products], [-1, 1]):
        for txt in string:
            identified_met = False
            lst = txt.strip(" ").split(" ")
            try:
                int(lst[0])
            except ValueError:
                lst = [" ".join(lst)]
            if (len(lst) == 2) and len(lst[0]):
                n_m = int(lst[0])
                m_name = lst[1].strip()
            else:
                n_m = 1
                m_name = lst[0].strip()

            if not len(m_name):
                continue
            
            try:
                m_ID = name_fix_dict_reverse[m_name] + "_c"
            except KeyError:
                m_ID = None

            for m in model.metabolites:
                if m.compartment == "c":
                    if m_ID:
                        if m.id == m_ID:
                            identified_met = True
                            m_dict[m] = k*n_m
                            break
                    else:
                        if m.name == m_name:
                            identified_met = True
                            m_dict[m] = k*n_m
                            break
            if not identified_met:
                m_dict[m_name] = k
                print("No met found for {0}, {1}".format(m_ID, m_name))
    return m_dict

def reaction_string_to_mets_name(string):
    substrate_string, product_string = string.split("->")
    substrates = substrate_string.split(" + ")
    products = product_string.split(" + ")

    m_dict = {}
    for string, k in zip([substrates, products], [-1, 1]):
        for txt in string:
            lst = txt.strip(" ").split(" ")
            try:
                int(lst[0])
            except ValueError:
                lst = [" ".join(lst)]
            if (len(lst) == 2) and len(lst[0]):
                n_m = int(lst[0])
                m_name = lst[1].strip()
            else:
                n_m = 1
                m_name = lst[0].strip()

            if not len(m_name):
                continue
            
            try:
                curr_value = m_dict[m_name]
            except KeyError:
                m_dict[m_name] = k*n_m
            else:
                m_dict[m_name] = curr_value + k*n_m
    return m_dict

def summarize_ahba_pathway(fn):
    df = pd.read_csv(fn, header = 0, sep = ";")
    m_sum_dict = {}
    for string in df["Reaction string"]:
        m_dict = reaction_string_to_mets_name(string)
        for m, value in m_dict.items():
            try:
                current_value = m_sum_dict[m]
            except KeyError:
                m_sum_dict[m] = value
            else:
                m_sum_dict[m] = current_value + value
    for key, value in m_sum_dict.items():
        if value != 0:
            print("{0}: {1}".format(key, value))


if __name__ == '__main__':
    if 0:
        all_pathways_to_table("../Data/constructed_pathways","../Data/constructed_pathways/csv/")
    if 0:
        fn = "../Data/validation_pathways"
        all_pathways_to_table(fn,fn+"/csv/")
    if 0:
        fix_BiGG_json()
    if 1:
        make_pathway_from_csv()
    if 0:
        fn = "../Data/ahba_synthesis.csv"
        summarize_ahba_pathway(fn)