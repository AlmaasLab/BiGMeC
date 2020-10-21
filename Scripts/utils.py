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
        
if __name__ == '__main__':
    all_pathways_to_table("../Data/constructed_pathways","../Data/constructed_pathways/csv/")
    # fn = r'C:\Users\snorres\OneDrive - SINTEF\Almaaslab\masteroppgaver\BGC-to-GEM\MasterThesis-master(1)\MasterThesis-master\gbk_db_output_models\gbk_db_output_models'
    # all_pathways_to_table(fn,fn+"/csv/")
    