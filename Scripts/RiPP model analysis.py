import os
from os.path import join
from pathlib import Path
from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model
import logging
from cobra import Model, Reaction, Metabolite
import cobra
from cobra.flux_analysis import gapfill
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion, flux_variability_analysis)
from cobra.medium import minimal_medium
import reframed
from cobra.util.solver import linear_reaction_coefficients
import cobra
from cobra.medium import minimal_medium
import pandas
import matplotlib.pyplot as plt
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import production_envelope
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

UNIMOD = cobra.io.read_sbml_model("Models\\Sco-GEM.xml")

MEDIUM3 = UNIMOD.medium.copy()

MEDIUM1 = {
        "EX_glc__D_e":10.0, "EX_nh4_e":1000.0, "EX_zn2_e":1000.0, "EX_pi_e":1000.0, "EX_h2o_e":1000.0, "EX_h_e":1000.0, "EX_mg2_e":1000.0, "EX_ca2_e":1000.0, "EX_o2_e":1000.0,
        "EX_k_e":1000.0, "EX_fe3_e":1000.0, "EX_fe2_e":1000.0, "EX_so4_e":1000.0, "EX_cl_e":1000.0, "EX_cobalt2_e":1000.0, "EX_cu2_e":1000.0, "EX_mn2_e":1000.0
    }

MEDIUM2 = {
    "EX_gly_e":10.0, "EX_ala__L_e":10.0, "EX_leu__L_e":10.0, "EX_met__L_e":10.0, "EX_phe__L_e":10.0, "EX_trp__L_e":10.0, "EX_lys__L_e":10.0, "EX_gln__L_e":10.0,
    "EX_glu__L_e":10.0, "EX_ser__L_e":10.0, "EX_pro__L_e":10.0, "EX_val__L_e":10.0, "EX_ile__L_e":10.0, "EX_cys__L_e":10.0, "EX_tyr__L_e":10.0, "EX_his__L_e":10.0,
    "EX_arg__L_e":10.0, "EX_asn__L_e":10.0, "EX_asp__L_e":10.0, "EX_thr__L_e":10.0, "EX_glc__D_e":10.0, "EX_nh4_e":1000.0, "EX_zn2_e":1000.0, "EX_pi_e":1000.0,
    "EX_h2o_e":1000.0, "EX_h_e":1000.0, "EX_mg2_e":1000.0, "EX_ca2_e":1000.0, "EX_o2_e":1000.0, "EX_k_e":1000.0, "EX_fe3_e":1000.0, "EX_fe2_e":1000.0, "EX_so4_e":1000.0,
    "EX_cl_e":1000.0, "EX_cobalt2_e":1000.0, "EX_cu2_e":1000.0, "EX_mn2_e":1000.0
}

TXID_ACCID_df = pd.read_csv(Path("Data\\TXID - GNMID.csv"))
BGCID_TXID_df = pd.read_csv(Path("Data\\RiPP BGCs - NCBI Tax ID.csv"))
GROWTH_FILTER_df = pd.read_csv("Data\\CarveMe extended model RiPP results.csv")

NCBI_ACCID_DICT = pd.Series(TXID_ACCID_df["Taxonomy ID"].values, index=TXID_ACCID_df["Assembly ID"]).to_dict()
TAXID_DICT = pd.Series(BGCID_TXID_df["BGC ID"].values, index=BGCID_TXID_df["Taxonomy ID"]).to_dict()

NCBI_ACCID_DICT_rev = pd.Series(TXID_ACCID_df["Assembly ID"].values, index=TXID_ACCID_df["Taxonomy ID"]).to_dict()
TAXID_DICT_rev = pd.Series(BGCID_TXID_df["Taxonomy ID"].values, index=BGCID_TXID_df["BGC ID"]).to_dict()

BGCID_GROWTH_FILTER_list = GROWTH_FILTER_df["BGC.ID"].values.tolist()
TAXID_GROWTH_FILTER_list = list(map(TAXID_DICT_rev.get, BGCID_GROWTH_FILTER_list))
ACCID_GROWTH_FILTER_list = list(map(NCBI_ACCID_DICT_rev.get, TAXID_GROWTH_FILTER_list))

def get_gb_list_from_antismash_output(cluster_path):  # yes
    # domains, pfam_entries, smCOG, EC_number, rxn_in, rxn_out, core_gene, amino_acid_sequence, predicted_EC, length
    # triple pound signs means that information is stored within a custom class object.
    gb_list = []
    for gb_record in SeqIO.parse(open(cluster_path, "r"), "genbank"):
        gb_list.append(gb_record)
    return gb_list

#Model analysis of either Sco-Gem or CarveMe models
def model_analysis(input):
    if input == "SRC":
        model_path = "RiPP Pathway Output\\CarveMe extended models"
    else:
        model_path = "RiPP Pathway Output\\Sco-GEM extended models"
    bgc_path = "Data\\MIBiG RiPP gbks"
    summary = pd.DataFrame({"BGC ID": [], "Core Number": [], "BGC Source Organism": [], "Core Class": [], "Core Subclass": [], "Max Growth": [], "Max NP": [], "Gradient": []})
    for file in tqdm(os.listdir(model_path)):
        new_row = {}
        model = read_sbml_model(os.path.join(model_path,file))
        med = MEDIUM3.copy()
        for k in MEDIUM3.keys():
            if k not in model.reactions.list_attr(attribute="id"):
                med.pop(k)
        model.medium = med
        O_id = "SK_"+os.path.splitext(file)[0].split("_")[1]+"_"+os.path.splitext(file)[0].split("_")[2]+"_mature_peptide_c"
        bgc_id = os.path.splitext(file)[0].split("_")[1]
        core_number = os.path.splitext(file)[0].split("_")[2]
        new_row["BGC ID"] = bgc_id
        new_row["Core Number"] = core_number
        new_row["Max Growth"] = model.slim_optimize()
        model.objective = O_id
        new_row["Max NP"] = model.slim_optimize()
        new_row["Gradient"] = None
        if new_row["Max Growth"] != 0:
            new_row["Gradient"] = new_row["Max NP"]/new_row["Max Growth"]
        gb_list = get_gb_list_from_antismash_output(os.path.join(bgc_path, bgc_id+".gbk"))
        new_row["BGC Source Organism"] = None
        new_row["Core Class"] = None
        new_row["Core Subclass"] = None
        for gb_record in gb_list:
            new_row["BGC Source Organism"] = gb_record.annotations["organism"]
            for feat in gb_record.features:
                if feat.type == "CDS_motif":
                    if feat.qualifiers.get("prepeptide") == ["core"]:
                        new_row["Core Class"] = feat.qualifiers.get("peptide")[0]
                        new_row["Core Subclass"] = feat.qualifiers.get("predicted_class")[0]
                        break
        temp_df = pd.DataFrame(new_row, index=[0])
        summary = pd.concat([temp_df, summary.loc[:]]).reset_index(drop=True)
    if input == "SRC":
        summary.to_csv("RiPP Pathway Output\\CarveMe model results.csv")

def heterologous_expression():
    model_path = "RiPP Pathway Output\\CarveMe heterologous extended models"
    bgc_path = "Data\\MIBiG RiPP gbks"
    summary = pd.DataFrame(
        {"BGC ID": [], "Core Number": [], "Host Model BGC": [], "Core Class": [], "Core Subclass": [], "Max Growth": [], "Max NP": [],
         "Gradient": []})
    for file in tqdm(os.listdir(model_path)):
        new_row = {}
        model = read_sbml_model(os.path.join(model_path, file))
        med = MEDIUM3.copy()
        for k in MEDIUM3.keys():
            if k not in model.reactions.list_attr(attribute="id"):
                med.pop(k)
        model.medium = med
        O_id = "DM_" + os.path.splitext(file)[0].split("_")[1] + "_" + os.path.splitext(file)[0].split("_")[2] + "_mature_peptide_c"
        bgc_id = os.path.splitext(file)[0].split("_")[1]
        core_number = os.path.splitext(file)[0].split("_")[2]
        new_row["BGC ID"] = bgc_id
        new_row["Core Number"] = core_number
        new_row["Max Growth"] = model.slim_optimize()
        model.objective = O_id
        new_row["Max NP"] = model.slim_optimize()
        new_row["Gradient"] = None
        if new_row["Max Growth"] != 0:
            new_row["Gradient"] = new_row["Max NP"] / new_row["Max Growth"]
        gb_list = get_gb_list_from_antismash_output(os.path.join(bgc_path, bgc_id + ".gbk"))
        new_row["Host Model BGC"] = TAXID_DICT[NCBI_ACCID_DICT[os.path.splitext(file)[0].split("_")[3]+"_"+os.path.splitext(file)[0].split("_")[4]]]
        new_row["Core Class"] = None
        new_row["Core Subclass"] = None
        for gb_record in gb_list:
            for feat in gb_record.features:
                if feat.type == "CDS_motif":
                    if feat.qualifiers.get("prepeptide") == ["core"]:
                        new_row["Core Class"] = feat.qualifiers.get("peptide")[0]
                        new_row["Core Subclass"] = feat.qualifiers.get("predicted_class")[0]
                        break
        temp_df = pd.DataFrame(new_row, index=[0])
        summary = pd.concat([temp_df, summary.loc[:]]).reset_index(drop=True)
    summary.to_csv("RiPP Pathway Output\\Heterologous model results.csv")

def main():
    model_analysis("SRC")
    #heterologous_expression()
if __name__ == "__main__":
    main()
