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

            substrate_string, product_string = r.reaction.split("-->")
            if not include_backbone:
                search_string = r"{0}+\_\d+".format(bgc_type)
                substrate_string = re.sub(r"\A\s\+|\+\s\Z", "", re.sub(search_string, "", substrate_string)).strip()
                product_string = re.sub(r"\A\s\+|\+\s\Z", "", re.sub(search_string, "", product_string)).strip()
                
            backbone_list.append([int(j), substrate_string, product_string])
    df = pd.DataFrame(backbone_list, columns = ["N", "Substrates", "Products"])
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

