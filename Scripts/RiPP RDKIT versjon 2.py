from cobra import Model, Reaction, Metabolite
from cobra.util.solver import linear_reaction_coefficients
import cobra
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion, flux_variability_analysis)
import os
from os.path import join
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import typing as t
from rdkit import Chem, DataStructs
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdChemReactions as Reactions
from rdkit.Chem.Draw import IPythonConsole
from PIL import Image
import matplotlib.pyplot as plt
import PIL
import pandas as pd
from Bio import SeqIO
from pathlib import Path

#Parameter for å avgjøre ringstørrelse
pRINGSIZE = 18
#Laster inn cobra GEM
#Dictionary for å gjøre peptidkode om til bigg stoff-id
biggdict = {"G":"gly_c", "A":"ala__L_c", "L":"leu__L_c", "M":"met__L_c", "F":"phe__L_c", "W":"trp__L_c", "K":"lys__L_c", "Q":"gln__L_c", "E":"glu__L_c", "S":"ser__L_c", "P":"pro__L_c", "V":"val__L_c", "I":"ile__L_c",
            "C":"cys__L_c", "Y":"tyr__L_c", "H":"his__L_c", "R":"arg__L_c", "N":"asn__L_c", "D":"asp__L_c", "T":"thr__L_c"}
#SMCOG1155 ansvarlig for dehydrering av serine/threonine ved glutamylering,
#LanB knyttet SMCOG1155, LanC tilknyttet SMCOG1140, LanM tilknyttet SMCOG1070, LanKC tilknyttet SMCOG1030, LanL tilknyttet SMCOG1030
#SMCOG1075 ansvarlig for kløyving av lederpeptid fra kjernepeptid
#SMCOG1030 ansvarlig for dehydrering av serine/threonine ved fosforylering + syklisering
#SMCOG1140 ansvarlig for syklisering mellom serine/threonine og cysteine
#LANC_like ansvarlig syklisering

#Node for å lage traverseringstre for reaksjonskjeder
class Metabolite_node:
    def __init__(self, metabolite, name, parent, reaction, molecule):
        self.metabolite = metabolite
        self.name = name
        self.parent = parent
        self.reaction = reaction
        self.molecule = molecule

#Funksjon for translasjon av peptidsekvens for alle RiPPs, returnerer initielle metabolittnode
def translation(peptidesequence, name, model):
    temp_mol = set_reactant_id(Chem.MolFromSequence(peptidesequence))
    Chem.SanitizeMol(temp_mol)
    smiles = Chem.MolToSmiles(temp_mol)
    prepeptide_reaction = Reaction(name+" prepeptide synthesis")
    prepeptide_reaction.name = "Peptide synthesis of "+name+" precursor"
    prepeptide_reaction.lower_bound = 0.
    prepeptide_reaction.upper_bound = 1000.

    prepeptide = Metabolite(
        id= name+"_prepeptide_c",
        formula= smiles,
        name=name+"precursor peptide",
        compartment="c"
    )

    aa_amount = 0
    for aa in peptidesequence:
        if aa in biggdict.keys():
            prepeptide_reaction.add_metabolites({model.metabolites.get_by_id(biggdict[aa]): -1.0})
            aa_amount += 1

    prepeptide_reaction.add_metabolites({model.metabolites.get_by_id("atp_c"):-aa_amount})
    prepeptide_reaction.add_metabolites({model.metabolites.get_by_id("amp_c"):aa_amount})
    prepeptide_reaction.add_metabolites({model.metabolites.get_by_id("ppi_c"):aa_amount})
    prepeptide_reaction.add_metabolites({model.metabolites.get_by_id("gtp_c"):-2*aa_amount})
    prepeptide_reaction.add_metabolites({model.metabolites.get_by_id("gdp_c"):2*aa_amount})
    prepeptide_reaction.add_metabolites({model.metabolites.get_by_id("pi_c"):2*aa_amount})
    prepeptide_reaction.add_metabolites({model.metabolites.get_by_id("h2o_c"):aa_amount-1})
    prepeptide_reaction.add_metabolites({prepeptide:1.0})

    m_node = Metabolite_node(prepeptide, name, None, prepeptide_reaction, temp_mol)

    return m_node

#Funksjon som kløyver kjerne- og ledepeptid fra hverandre
def core_cleave(node, coreseq, leaderseq, model):
    temp_mol1 = set_reactant_id(Chem.MolFromSequence(coreseq))
    Chem.SanitizeMol(temp_mol1)
    coresmiles = Chem.MolToSmiles(temp_mol1)
    temp_mol2 = set_reactant_id(Chem.MolFromSequence(leaderseq))
    Chem.SanitizeMol(temp_mol2)
    leadersmiles = Chem.MolToSmiles(temp_mol2)

    cleaved_leader_peptide = Metabolite(
        id=node.name+"_leader_peptide_c",
        formula=leadersmiles,
        name="Cleaved leader peptide of " + node.name,
        compartment="c"
    )

    cleaved_core_peptide = Metabolite(
        id=node.name+"_core_peptide_c",
        formula=coresmiles,
        name="Cleaved core peptide of " + node.name,
        compartment="c"
    )

    cleavage_reaction = Reaction("cleavage of " + node.name + " core and leader")
    cleavage_reaction.name = "Serine protease cleavage of " + node.name + " leader peptide from core peptide"
    cleavage_reaction.lower_bound = 0.
    cleavage_reaction.upper_bound = 1000.

    cleavage_reaction.add_metabolites(
        {node.metabolite: -1.0, model.metabolites.get_by_id("h2o_c"): -1.0, cleaved_core_peptide: 1.0,
         cleaved_leader_peptide: 1.0})

    m_node = Metabolite_node(cleaved_core_peptide, node.name, node, cleavage_reaction, temp_mol1)

    return m_node

#Må oppdateres til å bare ha med reagerbare områder (ubeskyttet)
def dehydratable(node):
    temp_mol = node.molecule
    serine_query = Chem.MolFromSmarts("[C:1]([CX4H2:2][OX2H:3])[C:4](=[O:5])")
    threonine_query = Chem.MolFromSmarts("[C:1]([C:2]([C:3])[O:4])[C:5](=[O:6])")
    #[C:1]([CX4H1:2]([CX4H3:3])[OX2H:4])[C:5](=[O:6])
    Chem.SanitizeMol(temp_mol)
    if temp_mol.HasSubstructMatch(serine_query) or temp_mol.HasSubstructMatch(threonine_query):
        return True
    else:
        return False

def cyclizationable(node):
    temp_mol = node.molecule
    serthr_query = Chem.MolFromSmarts("[C:6](=[C:7])[C:8](=[O:9])")
    cys_query = Chem.MolFromSmarts("[C:1]([C:2][SX2H:3])[C:4](=[O:5])")
    Chem.SanitizeMol(temp_mol)
    return temp_mol.HasSubstructMatch(serthr_query) and temp_mol.HasSubstructMatch(cys_query)

#Sjekker om det går an å gjøre N-term tailoring per LanO
def n_term_lac_able(node):
    temp_mol = node.molecule
    dha_n_term_query = Chem.MolFromSmarts("[NX3H2,NX2H1,NX4H3][C](=[C])[C](=[O])")
    Chem.SanitizeMol(temp_mol)
    return temp_mol.HasSubstructMatch(dha_n_term_query)

#NB bør testes i eget script. Funksjon som dehydrerer serine og threonine ved glutamylering
def glutamylation_elimination(node, model):
    serine_query = Chem.MolFromSmarts("[C:1]([CX4H2:2][OX2H:3])[C:4](=[O:5])")
    threonine_query = Chem.MolFromSmarts("[C:1]([C:2]([C:3])[O:4])[C:5](=[O:6])")
    substrate = node.molecule
    # Counts how many of each dehydration happens
    n = 0
    if substrate.HasSubstructMatch(serine_query):
        serine_smarts = "[C:1]([CX4H2:2][OX2H:3])[C:4](=[O:5]) >> [CX4H0:1](=[CX3H2:2])[C:4](=[O:5]).[O:3]"
        ser_rxn = AllChem.ReactionFromSmarts(serine_smarts)
        products = ser_rxn.RunReactants([substrate, ])

        while len(products) > 1:
            Chem.SanitizeMol(products[0][0])
            products = ser_rxn.RunReactants([products[0][0], ])
            Chem.SanitizeMol(products[0][0])
            n += 1

        substrate = products[0][0]

    if substrate.HasSubstructMatch(threonine_query):
        threonine_smarts = "[C:1]([C:2]([C:3])[O:4])[C:5](=[O:6]) >> [C:1](=[C:2][C:3])[C:5](=[O:6]).[O:4]"
        thr_rxn = AllChem.ReactionFromSmarts(threonine_smarts)
        products = thr_rxn.RunReactants([substrate, ])

        while len(products) > 1:
            Chem.SanitizeMol(products[0][0])
            products = thr_rxn.RunReactants([products[0][0], ])
            Chem.SanitizeMol(products[0][0])
            n += 1

        substrate = products[0][0]

    de_smiles = Chem.MolToSmiles(substrate)

    dehydrated_peptide = Metabolite(
        id=node.name+"_dehydrated_peptide_c",
        formula= de_smiles,
        name= node.name+" dehydrated peptide intermediate",
        compartment="c"
    )
    dehydration_reaction = Reaction("glut_elim_dehydration of " + node.name)
    dehydration_reaction.name = "Glutamylation elimination dehydration reaction of " + node.name
    dehydration_reaction.lower_bound = 0.
    dehydration_reaction.upper_bound = 1000.
    dehydration_reaction.add_metabolites({model.metabolites.get_by_id("glutrna_c"): -n, model.metabolites.get_by_id("glu__L_c"): n, model.metabolites.get_by_id("trnaglu_c"): n, node.metabolite: -1.0, dehydrated_peptide: 1.0})

    m_node = Metabolite_node(dehydrated_peptide, node.name, node, dehydration_reaction, products[0][0])
    return m_node

#NB bør testes i eget script. Funksjon som dehydrerer serine og threonine ved fosforylering
def phosphorylation_dehydration(node, model):
    serine_query = Chem.MolFromSmarts("[C:1]([CX4H2:2][OX2H:3])[C:4](=[O:5])")
    threonine_query = Chem.MolFromSmarts("[C:1]([C:2]([C:3])[O:4])[C:5](=[O:6])")
    substrate = node.molecule
    # Counts how many of each dehydration happens
    n = 0
    if substrate.HasSubstructMatch(serine_query):
        serine_smarts = "[C:1]([CX4H2:2][OX2H:3])[C:4](=[O:5]) >> [CX4H0:1](=[CX3H2:2])[C:4](=[O:5]).[O:3]"
        ser_rxn = AllChem.ReactionFromSmarts(serine_smarts)
        products = ser_rxn.RunReactants([substrate, ])

        while len(products) > 1:
            Chem.SanitizeMol(products[0][0])
            products = ser_rxn.RunReactants([products[0][0], ])
            Chem.SanitizeMol(products[0][0])
            n += 1

        substrate = products[0][0]

    if substrate.HasSubstructMatch(threonine_query):
        threonine_smarts = "[C:1]([C:2]([C:3])[O:4])[C:5](=[O:6]) >> [C:1](=[C:2][C:3])[C:5](=[O:6]).[O:4]"
        thr_rxn = AllChem.ReactionFromSmarts(threonine_smarts)
        products = thr_rxn.RunReactants([substrate, ])

        while len(products) > 1:
            Chem.SanitizeMol(products[0][0])
            products = thr_rxn.RunReactants([products[0][0], ])
            Chem.SanitizeMol(products[0][0])
            n += 1

        substrate = products[0][0]

    de_smiles = Chem.MolToSmiles(substrate)

    dehydrated_peptide = Metabolite(
        id=node.name + "_dehydrated_peptide_c",
        formula=de_smiles,
        name=node.name + " dehydrated peptide intermediate",
        compartment="c"
    )
    dehydration_reaction = Reaction("phospho_dehydration of " + node.name)
    dehydration_reaction.name = "Phosphorylation dehydration reaction of " + node.name
    dehydration_reaction.lower_bound = 0.
    dehydration_reaction.upper_bound = 1000.
    dehydration_reaction.add_metabolites(
        {model.metabolites.get_by_id("atp_c"): -n, model.metabolites.get_by_id("adp_c"): n,
         model.metabolites.get_by_id("pi_c"): n, node.metabolite: -1.0, dehydrated_peptide: 1.0})

    m_node = Metabolite_node(dehydrated_peptide, node.name, node, dehydration_reaction, products[0][0])
    return m_node

def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol

def mol_with_atom_rindex(mol):
    for atom in mol.GetAtoms():
        if "Reactant_idx" in atom.GetPropsAsDict():
                atom.SetAtomMapNum(atom.GetPropsAsDict()["Reactant_idx"])
    return mol

def set_reactant_id(mol):
    for atom in mol.GetAtoms():
        atom.SetIntProp('reactant_idx', atom.GetIdx())
        atom.SetAtomMapNum(atom.GetIdx())
    return mol

def map_reacted_residues(mol, queries):
    map_dict = {}
    for query in queries:
        m_query = Chem.MolFromSmarts(query)
        matches = mol.GetSubstructMatches(m_query)
        reacted = False
        count = -1
        direction_list = []
        for atom in mol.GetAtoms():
            reacted = False
            if any(atom.GetIdx() in l for l in matches):
                if "reactant_idx" in atom.GetPropsAsDict():
                    map_dict[atom.GetIdx()] = atom.GetPropsAsDict()["reactant_idx"]
                else:
                    map_dict[atom.GetIdx()] = None
                    reacted = True
                    count = 3
            if not reacted and count > 0 and "reactant_idx" in atom.GetPropsAsDict():
                direction_list.append(atom.GetPropsAsDict()["reactant_idx"])
                count -= 1
            if count == 0:
                direction = direction_list[0] - direction_list[1] + direction_list[1] - direction_list[2]
                if direction > 0:
                    pos = direction_list[0]+1
                    for ids in map_dict.keys():
                        if map_dict[ids] == None:
                            map_dict[ids] = pos
                            pos += 1
                elif direction < 0:
                    pos = direction_list[0]-1
                    for ids in map_dict.keys():
                        if map_dict[ids] == None:
                            map_dict[ids] = pos
                            pos -= 1
                else:
                    pos = direction_list[0]
                    for ids in map_dict.keys():
                        if map_dict[ids] == None:
                            map_dict[ids] = pos
                count -= 1
                direction_list = []
            elif count > 0 and atom.GetIdx() == len(mol.GetAtoms())-1:
                pos = direction_list[0]
                for ids in map_dict.keys():
                    if map_dict[ids] == None:
                        map_dict[ids] = pos
    return map_dict

def remap(matches, map):
    remapped = []
    for s in matches:
        site = []
        for p in s:
            site.append(map[p])
        remapped.append(site)
    return remapped

def average_position(list):
    new_list = []
    for l in list:
        new_list.append(sum(l)/len(l))
    new_list.sort()
    return new_list

def recursive_cyclization(mol, n):
    n += 1
    #SMARTS searching queries, first is for dehydrated ser/thr and second for cys
    queries = ["[C:6](=[C:7])[C:8](=[O:9])", "[C:1]([C:2][SX2H:3])[C:4](=[O:5])"]
    #Creates dictionary mapping current molecule idx to old value, with somewhat of the correct N-C direction conserved
    mol_map = map_reacted_residues(mol, queries)
    #Creates lists of reacting residues with average position using the remapped values
    mol_serthr = mol.GetSubstructMatches(Chem.MolFromSmarts(queries[0]))
    mol_cys = mol.GetSubstructMatches(Chem.MolFromSmarts(queries[1]))
    mol_cyslist = average_position(remap(mol_cys, mol_map))
    mol_serthrlist = average_position(remap(mol_serthr, mol_map))
    #Chooses residues best fitting the heuristic (cyclization in N-C direction, from cysteine to ser/thr. 1-2 residues minimum spacing for cycles. NB! 20 (18) is a temporary value, needs to check if this is reasonable)
    rightmost = mol_cyslist[len(mol_cyslist) - 1]
    optimal = mol_serthrlist[0]
    for st in mol_serthrlist:
        if abs(rightmost - st) < abs(rightmost - optimal) and st < rightmost - pRINGSIZE:
            optimal = st
    #Runs reaction
    smarts = "([C:1]([C:2][SX2H:3])[C:4](=[O:5]).[C:6](=[C:7])[C:8](=[O:9])) >> [C:1]([C:2][SX2H0:3][C:7][C:6][C:8](=[O:9]))[C:4](=[O:5])"
    cyc_rxn = AllChem.ReactionFromSmarts(smarts)
    products = cyc_rxn.RunReactants([mol, ])
    #Creates list containing remapped positional values of reacting cysteine residues for all of the products
    cycmol_cysmaps = []
    for p in products:
        Chem.SanitizeMol(p[0])
        cycmol_cys = p[0].GetSubstructMatches(Chem.MolFromSmarts(queries[1]))
        cycmol_map = map_reacted_residues(p[0], queries)
        cycmol_cysmaps.append(average_position(remap(cycmol_cys, cycmol_map)))
    #Creates list of products filtered by the cysteine condition (if cysteine is used up use remaining products instead)
    filtered_list = []
    if len(cycmol_cysmaps[0]) != 0:
        for m in range(len(cycmol_cysmaps)):
            if rightmost > cycmol_cysmaps[m][len(cycmol_cysmaps[m]) - 1]:
                filtered_list.append(products[m][0])
    else:
        for m in products:
            filtered_list.append(m[0])
    #Creates list containing remapped positional values for reacting ser/thr residues for the filtered products
    cycmol_serthrmaps = []
    for m in filtered_list:
        cycmol_serthr = m.GetSubstructMatches(Chem.MolFromSmarts(queries[0]))
        cycmol_map = map_reacted_residues(m, queries)
        cycmol_serthrmaps.append(average_position(remap(cycmol_serthr, cycmol_map)))
    #Finds index of the product best fitting the ser/thr condition
    i = 0
    dissim = 0
    for m in range(len(cycmol_serthrmaps)):
        temp_dissim = 0
        for st in cycmol_serthrmaps[m]:
            temp_dissim += abs(st - optimal)
        if temp_dissim > dissim:
            dissim = temp_dissim
            i = m

    cyc_mol = filtered_list[i]

    #image = Draw.MolToImage(cyc_mol, size=[2000, 2000])
    #image.show()
    if cyc_rxn.RunReactants([cyc_mol, ]) == ():
        return cyc_mol, n
    else:
        return recursive_cyclization(cyc_mol, n)

#Only supports lanthionine/methyllanthionine formation, labionin/avionin not implemented
def cyclization(node, model):
    substrate = node.molecule
    #Counts cycles (pun intended), each corresponds to creation of a water molecule
    n = 0
    product, c = recursive_cyclization(substrate, n)
    Chem.SanitizeMol(product)
    cyc_smiles = Chem.MolToSmiles(product)

    cyclic_peptide = Metabolite(
        id=node.name+"_cyclic_peptide_c",
        formula=cyc_smiles,
        name= node.name+" cyclic peptide",
        compartment="c"
    )

    cyclization_reaction = Reaction("cyclization of " + node.name)
    cyclization_reaction.name = "Cyclization reaction of " + node.name
    cyclization_reaction.lower_bound = 0.
    cyclization_reaction.upper_bound = 1000.
    cyclization_reaction.add_metabolites(
        {node.metabolite:-1.0, cyclic_peptide:1.0, model.metabolites.get_by_id("h2o_c"):c}
    )

    m_node = Metabolite_node(cyclic_peptide, node.name, node, cyclization_reaction, product)

    return m_node

def lan_o_tailoring(node, model):
    substrate = node.molecule
    dha_lac_reduction_smarts = "[NX3H2,NX2H1,NX4H3:1][C:2](=[C:3])[C:4](=[O:5]) >> [OH][C:2]([C:3])[C:4](=[O:5]).[N:1]"
    dlr_rxn = AllChem.ReactionFromSmarts(dha_lac_reduction_smarts)
    product = dlr_rxn.RunReactants([substrate,])

    Chem.SanitizeMol(product[0][0])
    substrate = product[0][0]

    lan_o_smiles = Chem.MolToSmiles(substrate)

    lan_o_tailored_peptide = Metabolite(
        id=node.name + "_lan_o_tailored_peptide_c",
        formula=lan_o_smiles,
        name=node.name + " LanO tailored peptide",
        compartment="c"
    )

    lan_o_tailoring_reaction = Reaction("lan_o_tailoring_reaction of "+node.name)
    lan_o_tailoring_reaction.name = "LanO tailoring reaction of "+node.name
    lan_o_tailoring_reaction.lower_bound = 0.
    lan_o_tailoring_reaction.upper_bound = 1000.
    lan_o_tailoring_reaction.add_metabolites(
        {node.metabolite:-1.0, lan_o_tailored_peptide:1.0, model.metabolites.get_by_id("h2o_c"):-1.0, model.metabolites.get_by_id("nadph_c"):-1.0,
         model.metabolites.get_by_id("nadp_c"):1.0, model.metabolites.get_by_id("nh4_c"):1.0}
    )

    m_node = Metabolite_node(lan_o_tailored_peptide, node.name, node, lan_o_tailoring_reaction, substrate)
    return m_node

def get_gb_list_from_antismash_output(cluster_path):  # yes
    # domains, pfam_entries, smCOG, EC_number, rxn_in, rxn_out, core_gene, amino_acid_sequence, predicted_EC, length
    # triple pound signs means that information is stored within a custom class object.
    gb_list = []
    for gb_record in SeqIO.parse(open(cluster_path, "r"), "genbank"):
        gb_list.append(gb_record)
    return gb_list

def run(bgc_path, output):
    print("Running...")
    summary = pd.DataFrame({"Name":[], "BGC":[], "Success":[], "End product SMILES":[], "Metabolite list":[]})
    column_list = ["Name", "BGC", "Success", "End product SMILES", "Metabolite list"]
    bgc_path = Path(bgc_path)

    if bgc_path.is_dir():
        for filename in os.listdir(bgc_path):
            p = os.path.join(bgc_path, filename)
            if os.path.isfile(p):
                print(filename)
                bgc_dict = {"Name":os.path.splitext(filename)[0], "BGC":"Unknown", "Success":"Failed", "End product SMILES":"Unknown", "Metabolite list":[],
                        "SMCOG1070":False, "SMCOG1140":False, "SMCOG1155":False, "SMCOG1075":False, "SMCOG1030":False, "Core":None, "Leader":None, "LanO":False}
                result = _run(p, output, bgc_dict)
                new_row = {}
                for column in column_list:
                    new_row[column] = [result[column]]
                temp_df = pd.DataFrame(new_row)
                summary = pd.concat([temp_df, summary.loc[:]]).reset_index(drop=True)
    else:
        name = os.path.basename(bgc_path)
        bgc_dict = {"Name": os.path.splitext(name)[0], "BGC": "Unknown", "Success": "Failed", "End product SMILES": "Unknown", "Metabolite list": [],
                    "SMCOG1070": False, "SMCOG1140": False, "SMCOG1155": False, "SMCOG1075": False, "SMCOG1030":False, "Core": None, "Leader": None}
        result = _run(bgc_path, output, bgc_dict)
        new_row = {}
        for column in column_list:
            new_row[column] = [result[column]]
        temp_df = pd.DataFrame(new_row)
        summary = pd.concat([temp_df, summary.loc[:]]).reset_index(drop=True)

    summary.to_csv(os.path.join(output, "test_med_tailoring.csv"))

def parse_gbk(bgc_path, bgc_dict):
    gb_list = get_gb_list_from_antismash_output(bgc_path)
    for gb_record in gb_list:
        for feat in gb_record.features:
            if feat.type == "CDS_motif":
                if feat.qualifiers.get("prepeptide") == ["core"]:
                    bgc_dict["Core"] = feat.qualifiers.get("core_sequence")
                    bgc_dict["Leader"] = feat.qualifiers.get("leader_sequence")
            elif feat.type == "CDS":
                if feat.qualifiers.get("gene_functions") != None:
                    for i in feat.qualifiers.get("gene_functions"):
                        if "smcogs" in i:
                            pos = i.find("smcogs") + 8
                            bgc_dict[i[pos:pos + 9]] = True
                            if feat.qualifiers.get("gene") != None:
                                if i[pos:pos + 9] in ["SMCOG1007", "SMCOG1001", "SMCOG1251"] and feat.qualifiers.get("gene")[0][3:4] == "O":
                                    bgc_dict["LanO"] = True
            elif feat.type == "misc_feature":
                bgc_dict["BGC"] = feat.qualifiers.get("note")
    return bgc_dict

def _run(bgc_path, output, bgc_dict):
    bgc_dict = parse_gbk(bgc_path, bgc_dict)
    if bgc_dict["Core"] == None and bgc_dict["Leader"] == None:
        bgc_dict["Success"] = "Missing core and leader peptide sequence"
        return bgc_dict
    elif bgc_dict["Core"] == None:
        bgc_dict["Success"] = "Missing core peptide sequence"
        return bgc_dict
    elif bgc_dict["Leader"] == None:
        bgc_dict["Success"] = "Missing leader peptide sequence"
        return bgc_dict
    #Run main script here using altered dictionary, main script saves sbml and PIL files in same output folder
    bgc_dict = add_to_model(bgc_dict, output)
    return bgc_dict

def add_to_model(bgc_dict, output):
    model = cobra.io.read_sbml_model(
        join("C:\\Users\\adtho\\Documents\\Python Scripts\\pythonProject-20221219T143334Z-001\\pythonProject\\models",
             "SCO-GEM.xml"))

    leaderseq = bgc_dict["Leader"][0]
    coreseq = bgc_dict["Core"][0]
    pepseq = leaderseq+coreseq

    origin_node = translation(pepseq, bgc_dict["Name"], model)
    current_node = origin_node
    current_node = core_cleave(current_node, coreseq, leaderseq, model)
    if dehydratable(current_node) and bgc_dict["SMCOG1155"]:
        current_node = glutamylation_elimination(current_node, model)
    if dehydratable(current_node) and bgc_dict["SMCOG1030"] or dehydratable(current_node) and bgc_dict["SMCOG1070"]:
        current_node = phosphorylation_dehydration(current_node, model)
    if cyclizationable(current_node) and bgc_dict["SMCOG1140"] or cyclizationable(current_node) and bgc_dict["SMCOG1070"] or cyclizationable(current_node) and bgc_dict["SMCOG1030"]:
        current_node = cyclization(current_node, model)
    if n_term_lac_able(current_node) and bgc_dict["LanO"]:
        current_node = lan_o_tailoring(current_node, model)

    bgc_dict["End product SMILES"] = current_node.metabolite.formula
    metabolite_list = []

    while current_node != None:
        #print(current_node.reaction)
        metabolite_list.append(current_node.metabolite.id)
        current_node = current_node.parent
    metabolite_list.reverse()

    bgc_dict["Metabolite list"] = metabolite_list
    bgc_dict["Success"] = "Yes"
    #Missing functionality: Output sbml file for model, some sort of image handling with appropriate naming, fingerprint molecular matching for characterized products
    return bgc_dict

def main():
    run("C:\\Users\\adtho\\Documents\\Python Scripts\\pythonProject-20221219T143334Z-001\\pythonProject\\antiSMASH 6 annotated RiPP clusters (lanthipeptides)",
        "C:\\Users\\adtho\\Documents\\Python Scripts\\pythonProject-20221219T143334Z-001\\pythonProject\\Output")

if __name__ == "__main__":
    main()