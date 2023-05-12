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
from rdkit.Chem import Draw, rdCoordGen
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdChemReactions as Reactions
from rdkit.Chem.Draw import IPythonConsole
from PIL import Image
import matplotlib.pyplot as plt
import PIL
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from tqdm import tqdm
os.chdir('..')

TXID_ACCID_df = pd.read_csv(Path("Data\\TXID - GNMID.csv"))
BGCID_TXID_df = pd.read_csv(Path("Data\\RiPP BGCs - NCBI Tax ID.csv"))
BGCID_SMILES_df = pd.read_csv(Path("Data\\BGCID - PubChem SMILES.csv"))

NCBI_ACCID_DICT = pd.Series(TXID_ACCID_df["Assembly ID"].values, index=TXID_ACCID_df["Taxonomy ID"]).to_dict()
TAXID_DICT = pd.Series(BGCID_TXID_df["Taxonomy ID"].values, index=BGCID_TXID_df["BGC ID"]).to_dict()
BGC_to_smiles = pd.Series(BGCID_SMILES_df["SMILES"].values, index=BGCID_SMILES_df["BGC ID"]).to_dict()

UNIMOD = cobra.io.read_sbml_model("Models\\Sco-GEM.xml")

#Parameter deciding ring size in thioether crosslinking reactions (lanthipeptides mostly)
pRINGSIZE = 9
#Dictionary of BIGG compound IDs for each amino acid code
biggdict = {"G":"gly_c", "A":"ala__L_c", "L":"leu__L_c", "M":"met__L_c", "F":"phe__L_c", "W":"trp__L_c", "K":"lys__L_c", "Q":"gln__L_c", "E":"glu__L_c", "S":"ser__L_c", "P":"pro__L_c", "V":"val__L_c", "I":"ile__L_c",
            "C":"cys__L_c", "Y":"tyr__L_c", "H":"his__L_c", "R":"arg__L_c", "N":"asn__L_c", "D":"asp__L_c", "T":"thr__L_c"}
#SMCOG1155 ansvarlig for dehydrering av serine/threonine ved glutamylering,
#LanB knyttet SMCOG1155, LanC tilknyttet SMCOG1140, LanM tilknyttet SMCOG1070, LanKC tilknyttet SMCOG1030, LanL tilknyttet SMCOG1030
#SMCOG1075 ansvarlig for kløyving av lederpeptid fra kjernepeptid
#SMCOG1030 ansvarlig for dehydrering av serine/threonine ved fosforylering + syklisering
#SMCOG1140 ansvarlig for syklisering mellom serine/threonine og cysteine
#LANC_like ansvarlig syklisering

#This class makes node objects which represent each metabolite in the reaction pathway
#Class of objects for associating COBRApy metabolite and reaction objects as well as RDKit molecule object, SMILES string and metabolite name to a given metabolite in the RiPP pathway
class Metabolite_node:
    def __init__(self, metabolite, name, parent, reaction, molecule, smiles):
        self.metabolite = metabolite
        self.name = name
        self.parent = parent
        self.reaction = reaction
        self.molecule = molecule
        self.smiles = smiles

#Function for the translation reaction of the precursor peptide, returns the initial metabolite node
def translation(peptidesequence, coresequence, name, model):
    temp_mol = set_reactant_id(Chem.MolFromSequence(coresequence))
    Chem.SanitizeMol(temp_mol)
    smiles = Chem.MolToSmiles(temp_mol)
    prepeptide_reaction = Reaction(name+"_prepeptide_synthesis")
    prepeptide_reaction.name = "Peptide synthesis of "+name+" precursor"
    prepeptide_reaction.lower_bound = 0.
    prepeptide_reaction.upper_bound = 1000.

    prepeptide = Metabolite(
        id= name+"_prepeptide_c",
        formula= AllChem.CalcMolFormula(temp_mol),
        name=name+"precursor peptide",
        compartment="c"
    )

    aa_amount = 0
    for aa in peptidesequence:
        if aa in biggdict.keys():
            prepeptide_reaction.add_metabolites({UNIMOD.metabolites.get_by_id(biggdict[aa]): -1.0})
            aa_amount += 1

    prepeptide_reaction.add_metabolites({UNIMOD.metabolites.get_by_id("atp_c"):-aa_amount})
    prepeptide_reaction.add_metabolites({UNIMOD.metabolites.get_by_id("amp_c"):aa_amount})
    prepeptide_reaction.add_metabolites({UNIMOD.metabolites.get_by_id("ppi_c"):aa_amount})
    prepeptide_reaction.add_metabolites({UNIMOD.metabolites.get_by_id("gtp_c"):-2*aa_amount})
    prepeptide_reaction.add_metabolites({UNIMOD.metabolites.get_by_id("gdp_c"):2*aa_amount})
    prepeptide_reaction.add_metabolites({UNIMOD.metabolites.get_by_id("pi_c"):2*aa_amount})
    prepeptide_reaction.add_metabolites({UNIMOD.metabolites.get_by_id("h2o_c"):aa_amount-1})
    prepeptide_reaction.add_metabolites({prepeptide:1.0})

    m_node = Metabolite_node(prepeptide, name, None, prepeptide_reaction, temp_mol, smiles)
    return m_node

#Function for the peptides cleavage reaction of core and leader peptides in lanthi/lassopeptides
def core_cleave(node, leaderseq, model, bgc_dict):
    temp_mol1 = node.molecule
    Chem.SanitizeMol(temp_mol1)
    coresmiles = Chem.MolToSmiles(temp_mol1)
    temp_mol2 = set_reactant_id(Chem.MolFromSequence(leaderseq))
    Chem.SanitizeMol(temp_mol2)
    leadersmiles = Chem.MolToSmiles(temp_mol2)

    cleaved_leader_peptide = Metabolite(
        id=node.name+"_leader_peptide_c",
        formula=AllChem.CalcMolFormula(temp_mol2),
        name="Cleaved leader peptide of " + node.name,
        compartment="c"
    )

    cleaved_core_peptide = Metabolite(
        id=node.name+"_core_peptide_c",
        formula=AllChem.CalcMolFormula(temp_mol1),
        name="Cleaved core peptide of " + node.name,
        compartment="c"
    )

    cleavage_reaction = Reaction("cleavage_of_" + node.name + "_core_and_leader")
    cleavage_reaction.name = "Peptidase cleavage of " + node.name + " leader peptide from core peptide"
    cleavage_reaction.lower_bound = 0.
    cleavage_reaction.upper_bound = 1000.
    if bgc_dict["Core Class"] == "lanthipeptide":
        cleavage_reaction.add_metabolites(
            {node.metabolite: -1.0, UNIMOD.metabolites.get_by_id("h2o_c"): -1.0, cleaved_core_peptide: 1.0})
    if bgc_dict["Core Class"] == "lassopeptide":
        cleavage_reaction.add_metabolites(
            {node.metabolite: -1.0, UNIMOD.metabolites.get_by_id("h2o_c"): -1.0, UNIMOD.metabolites.get_by_id("atp_c"):-1.0,cleaved_core_peptide: 1.0,
             UNIMOD.metabolites.get_by_id("adp_c"):1.0, UNIMOD.metabolites.get_by_id("pi_c"):1.0})
    m_node = Metabolite_node(cleaved_core_peptide, node.name, node, cleavage_reaction, temp_mol1, coresmiles)

    return m_node

def transport(node, model):

    mature_secreted_peptide = Metabolite(
        id=node.name + "_mature_peptide_e",
        formula=node.metabolite.formula,
        name="Mature peptide of " + node.name,
        compartment="e"
    )

    transport_reaction = Reaction("transport_of_mature_"+node.name+"_peptide")
    transport_reaction.name = "Transport of "+node.name+" mature peptide"
    transport_reaction.lower_bound = 0.
    transport_reaction.upper_bound = 1000.

    transport_reaction.add_metabolites(
        {node.metabolite:-1.0, mature_secreted_peptide:1.0}
    )

    m_node = Metabolite_node(mature_secreted_peptide, node.name, node, transport_reaction, node.molecule, node.formula)

    return m_node

#Checks if metabolite has residues (Ser/Thr) which can be dehydrated by associated enzymes
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

#Checks if metabolite contains chemical motifs (Dhx + Cys) which can partake in thioether crosslinking
def cyclizationable(node):
    temp_mol = node.molecule
    serthr_query = Chem.MolFromSmarts("[C:6](=[C:7])[C:8](=[O:9])")
    cys_query = Chem.MolFromSmarts("[C:1]([C:2][SX2H:3])[C:4](=[O:5])")
    Chem.SanitizeMol(temp_mol)
    return temp_mol.HasSubstructMatch(serthr_query) and temp_mol.HasSubstructMatch(cys_query)

#Checks if it is possible to do N-terminus tailoring per LanO activity
def n_term_lac_able(node):
    temp_mol = node.molecule
    dha_n_term_query = Chem.MolFromSmarts("[NX3H2,NX2H1,NX4H3][C](=[C])[C](=[O])")
    Chem.SanitizeMol(temp_mol)
    return temp_mol.HasSubstructMatch(dha_n_term_query)

#Checks if it is possible to do C-terminus tailoring per LanD activity
def c_term_oxidative_decarboxylationable(node):
    temp_mol = node.molecule
    cys_c_term_query = Chem.MolFromSmarts("[NX3H1][C]([C][S])[C](=[O])[OX2H1]")
    Chem.SanitizeMol(temp_mol)
    return temp_mol.HasSubstructMatch(cys_c_term_query)

#Checks if it is possible to perform heterocyclization of cysteine residues (thiopeptides)
def heterocyclizationable(node):
    temp_mol = node.molecule
    cys_query = Chem.MolFromSmarts("[C:1](=[O:2])[N:3][C:4]([C:5][SX2H:6])[C:7](=[O:8])")
    Chem.SanitizeMol(temp_mol)
    return temp_mol.HasSubstructMatch(cys_query)

#Checks if metabolite contains chemical motifs (2-pi + 4-pi) necessary for thiopeptide macrocyclization
def macrocyclizationable(node):
    temp_mol = node.molecule
    two_pi_query = Chem.MolFromSmarts("[N:1][C:2](=[CX3H2:3])")
    four_pi_query = Chem.MolFromSmarts("[C:9](=[CX3H2:10])([c:11]2[s:12][c:13][c:14][n:15]2)[N:16][C:17](=[O:18])[c:19]3[n:20][c:21][s:22][c:23]3")
    Chem.SanitizeMol(temp_mol)
    return temp_mol.HasSubstructMatch(two_pi_query) and temp_mol.HasSubstructMatch(four_pi_query)

#Checks if metabolite contains chemical motifs for performing macrolactam ring cyclization
def macrolactamable(node):
    temp_mol = node.molecule
    carboxy_query = Chem.MolFromSmarts("[CX4H2][CX3](=[O])[OX2H]")
    n_term_query = Chem.MolFromSmarts("[NX3H2][CX4]")
    Chem.SanitizeMol(temp_mol)
    return temp_mol.HasSubstructMatch(carboxy_query) and temp_mol.HasSubstructMatch(n_term_query)

#Checks if metabolite contains more than 1 cysteine residue
def disulfide_bridgable(node):
    temp_mol = node.molecule
    cys_query = Chem.MolFromSmarts("[C:1]([C:2][SX2H:3])[C:4](=[O:5])")
    Chem.SanitizeMol(temp_mol)
    return len(temp_mol.GetSubstructMatches(cys_query)) > 1

#Checks if it is possible to perform C-terminus thioenol to thioether crosslinking
def contains_thioenol_c_term(mol):
    temp_mol = mol
    thioenol_c_term_query = Chem.MolFromSmarts("[NX3H1:1][C:2]=[C:3][SX2H1:4]")
    Chem.SanitizeMol(temp_mol)
    return temp_mol.HasSubstructMatch(thioenol_c_term_query)

#Function that dehydrates serine and threonine per the glutamylation elimination reaction catalyzed by LanB (class I lanthipeptides)
def glutamylation_elimination(node, model, strictness):
    if strictness == 0:
        serine_query = Chem.MolFromSmarts("[C:1]([CX4H2:2][OX2H:3])[C:4](=[O:5])")
        threonine_query = Chem.MolFromSmarts("[C:1]([C:2]([C:3])[O:4])[C:5](=[O:6])")
        serine_smarts = "[C:1]([CX4H2:2][OX2H:3])[C:4](=[O:5]) >> [CX4H0:1](=[CX3H2:2])[C:4](=[O:5]).[O:3]"
        threonine_smarts = "[C:1]([C:2]([C:3])[O:4])[C:5](=[O:6]) >> [C:1](=[C:2][C:3])[C:5](=[O:6]).[O:4]"
    elif strictness == 1:
        serine_query = Chem.MolFromSmarts("[C:1]([CX4H2:2][OX2H:3])[C,c:4]")
        threonine_query = Chem.MolFromSmarts("[C:1]([C:2]([C:3])[OX2H:4])[C,c:5]")
        serine_smarts = "[C:1]([CX4H2:2][OX2H:3])[C,c:4] >> [CX4H0:1](=[CX3H2:2])[C,c:4].[O:3]"
        threonine_smarts = "[C:1]([C:2]([C:3])[OX2H:4])[C,c:5] >> [C:1](=[C:2][C:3])[C,c:5].[O:4]"
    else:
        serine_query = Chem.MolFromSmarts("[NX3H2,NX3H1,NX3H0,NX2H1,NX4H3:1][CX4H1:2][CX4H2,CX4H1:3][O:4]")
        threonine_query = Chem.MolFromSmarts("[NX3H2,NX3H1,NX3H0,NX2H1,NX4H3:1][CX4H1:2][CX4H2,CX4H1:3][O:4]")
        serine_smarts = "[NX3H2,NX3H1,NX3H0:1][CX4H1:2][CX4H2,CX4H1:3][O:4] >> [NX3:1][CX3H0:2]=[CX3:3].[O:4]"
        threonine_smarts = "[NX3H2,NX3H1,NX3H0:1][CX4H1:2][CX4H2,CX4H1:3][O:4] >> [NX3:1][CX3H0:2]=[CX3:3].[O:4]"

    substrate = node.molecule
    # Counts how many of each dehydration happens
    n = 0
    if substrate.HasSubstructMatch(serine_query):
        ser_rxn = AllChem.ReactionFromSmarts(serine_smarts)
        products = ser_rxn.RunReactants([substrate, ])
        n += 1
        while len(products) > 1:
            Chem.SanitizeMol(products[0][0])
            products = ser_rxn.RunReactants([products[0][0], ])
            Chem.SanitizeMol(products[0][0])
            n += 1

        substrate = products[0][0]

    if substrate.HasSubstructMatch(threonine_query):
        thr_rxn = AllChem.ReactionFromSmarts(threonine_smarts)
        products = thr_rxn.RunReactants([substrate, ])
        n += 1
        while len(products) > 1:
            Chem.SanitizeMol(products[0][0])
            products = thr_rxn.RunReactants([products[0][0], ])
            Chem.SanitizeMol(products[0][0])
            n += 1

        substrate = products[0][0]

    AllChem.SanitizeMol(substrate)
    de_smiles = Chem.MolToSmiles(substrate)

    dehydrated_peptide = Metabolite(
        id=node.name+"_dehydrated_peptide_c",
        formula= AllChem.CalcMolFormula(substrate),
        name= node.name+" dehydrated peptide intermediate",
        compartment="c"
    )
    dehydration_reaction = Reaction("glut_elim_dehydration_of_" + node.name)
    dehydration_reaction.name = "Glutamylation elimination dehydration reaction of " + node.name
    dehydration_reaction.lower_bound = 0.
    dehydration_reaction.upper_bound = 1000.
    dehydration_reaction.add_metabolites({UNIMOD.metabolites.get_by_id("glutrna_c"): -n, UNIMOD.metabolites.get_by_id("glu__L_c"): n, UNIMOD.metabolites.get_by_id("trnaglu_c"): n, node.metabolite: -1.0, dehydrated_peptide: 1.0})

    m_node = Metabolite_node(dehydrated_peptide, node.name, node, dehydration_reaction, substrate, de_smiles)
    return m_node

#Function that dehydrates serine and threonine per the phosphorylation reaction catalyzed by LanM/KC/L (class II-IV lanthipeptides)
def phosphorylation_dehydration(node, model, strictness):
    if strictness == 0:
        serine_query = Chem.MolFromSmarts("[C:1]([CX4H2:2][OX2H:3])[C:4](=[O:5])")
        threonine_query = Chem.MolFromSmarts("[C:1]([C:2]([C:3])[O:4])[C:5](=[O:6])")
        serine_smarts = "[C:1]([CX4H2:2][OX2H:3])[C:4](=[O:5]) >> [CX4H0:1](=[CX3H2:2])[C:4](=[O:5]).[O:3]"
        threonine_smarts = "[C:1]([C:2]([C:3])[O:4])[C:5](=[O:6]) >> [C:1](=[C:2][C:3])[C:5](=[O:6]).[O:4]"
    elif strictness == 1:
        serine_query = Chem.MolFromSmarts("[C:1]([C:2][O:3])[C:4](=[O:5])")
        threonine_query = Chem.MolFromSmarts("[C:1]([C:2]([C:3])[O:4])[C:5](=[O:6])")
        serine_smarts = "[C:1]([C:2][O:3])[C:4](=[O:5]) >> [C:1](=[C:2])[C:4](=[O:5]).[O:3]"
        threonine_smarts = "[C:1]([C:2]([C:3])[O:4])[C:5](=[O:6]) >> [C:1](=[C:2][C:3])[C:5](=[O:6]).[O:4]"
    else:
        serine_query = Chem.MolFromSmarts("[C:1]([C:2][O:3])[C:4](=[O:5])")
        threonine_query = Chem.MolFromSmarts("[C:1]([C:2]([C:3])[O:4])[C:5](=[O:6])")
        serine_smarts = "[C:1]([C:2][O:3])[C:4](=[O:5]) >> [C:1](=[C:2])[C:4](=[O:5]).[O:3]"
        threonine_smarts = "[C:1]([C:2]([C:3])[O:4])[C:5](=[O:6]) >> [C:1](=[C:2][C:3])[C:5](=[O:6]).[O:4]"
    substrate = node.molecule
    # Counts how many of each dehydration happens
    n = 0
    if substrate.HasSubstructMatch(serine_query):
        ser_rxn = AllChem.ReactionFromSmarts(serine_smarts)
        products = ser_rxn.RunReactants([substrate, ])
        n += 1
        while len(products) > 1:
            Chem.SanitizeMol(products[0][0])
            products = ser_rxn.RunReactants([products[0][0], ])
            Chem.SanitizeMol(products[0][0])
            n += 1

        substrate = products[0][0]
    if substrate.HasSubstructMatch(threonine_query):
        thr_rxn = AllChem.ReactionFromSmarts(threonine_smarts)
        products = thr_rxn.RunReactants([substrate, ])
        n += 1
        while len(products) > 1:
            Chem.SanitizeMol(products[0][0])
            products = thr_rxn.RunReactants([products[0][0], ])
            Chem.SanitizeMol(products[0][0])
            n += 1

        substrate = products[0][0]
    AllChem.SanitizeMol(substrate)
    de_smiles = Chem.MolToSmiles(substrate)

    dehydrated_peptide = Metabolite(
        id=node.name + "_dehydrated_peptide_c",
        formula=AllChem.CalcMolFormula(substrate),
        name=node.name + " dehydrated peptide intermediate",
        compartment="c"
    )
    dehydration_reaction = Reaction("phospho_dehydration_of_" + node.name)
    dehydration_reaction.name = "Phosphorylation dehydration reaction of " + node.name
    dehydration_reaction.lower_bound = 0.
    dehydration_reaction.upper_bound = 1000.
    dehydration_reaction.add_metabolites(
        {UNIMOD.metabolites.get_by_id("atp_c"): -n, UNIMOD.metabolites.get_by_id("adp_c"): n,
         UNIMOD.metabolites.get_by_id("pi_c"): n, node.metabolite: -1.0, dehydrated_peptide: 1.0})

    m_node = Metabolite_node(dehydrated_peptide, node.name, node, dehydration_reaction, products[0][0], de_smiles)
    return m_node

#Sets up the reactant ids for the atoms of a newly created polypeptide, so that the id represents the original N-C direction
def set_reactant_id(mol):
    for atom in mol.GetAtoms():
        atom.SetIntProp('reactant_idx', atom.GetIdx())
    return mol

#Maps the id of atoms which have participated in reactions (and therefore lost their reactant ids), using inference by the presence of reactant ids in nearby reactions
def map_reacted_residues(mol, queries):
    map_dict = {}
    for query in queries:
        m_query = Chem.MolFromSmarts(query)
        matches = mol.GetSubstructMatches(m_query)
        count = -1
        direction_list = []
        emergency_pos = 0
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
                    emergency_pos = pos
                    for ids in map_dict.keys():
                        if map_dict[ids] == None:
                            map_dict[ids] = pos
                            pos += 1
                elif direction < 0:
                    pos = direction_list[0]-1
                    emergency_pos = pos
                    for ids in map_dict.keys():
                        if map_dict[ids] == None:
                            map_dict[ids] = pos
                            pos -= 1
                else:
                    pos = direction_list[0]
                    emergency_pos = pos
                    for ids in map_dict.keys():
                        if map_dict[ids] == None:
                            map_dict[ids] = pos
                count -= 1
                direction_list = []
            elif count > 0 and atom.GetIdx() == len(mol.GetAtoms())-1:
                if len(direction_list) == 0 and "reactant_idx" in atom.GetPropsAsDict():
                    direction_list.append(atom.GetPropsAsDict()["reactant_idx"])
                if len(direction_list) == 0:
                    #Best estimate in case there is no reactant idx atoms near the unmapped atoms for a while, probably inaccurate
                    direction_list.append(emergency_pos)
                pos = direction_list[0]
                for ids in map_dict.keys():
                    if map_dict[ids] == None:
                        map_dict[ids] = pos
    return map_dict

#Remaps the position of residues found through substructure searching (which returns atom_id), using the mapped reactant ids
def remap(matches, map):
    remapped = []
    for s in matches:
        site = []
        for p in s:
            site.append(map[p])
        remapped.append(site)
    return remapped

#Calculates an average position value for a substructure
def average_position(list):
    new_list = []
    for l in list:
        new_list.append(sum(l)/len(l))
    new_list.sort()
    return new_list

#Models thioether crosslinking between cysteine/thioenolate and Dhx, returns optimal fully cyclizied product
def recursive_cyclization(mol, n, bgc_dict, mode):
    #Sets up standard searching queries/parameters (for regular lanthipeptide cyclization)
    cysquery = 1
    dhxquery = 0
    ringsize = pRINGSIZE
    #Counts recursive cycles
    n += 1
    # SMARTS searching queries, first is for dehydrated ser/thr, second for cys, and third for c-term thioenolate
    queries = ["[C:6](=[C:7])[C:8](=[O:9])", "[C:1]([C:2][SX2H:3])[C:4](=[O:5])", "[NX3H1:1][C:2]=[C:3][SX2H1:4]", "[C:6](=[CX3H2:7])[N:8][C:9](=[O:10])", "[C:6](=[C:7])[N:8][C:9](=[O:10])"]
    # Creates dictionary mapping current molecule idx to old value, with somewhat of the correct N-C direction conserved
    mol_map = map_reacted_residues(mol, queries)
    #Chooses reaction between thioenol attack or cysteine attack (thioenol intermediate is a C-term tailoring done by LanD)
    if contains_thioenol_c_term(mol) and mol.HasSubstructMatch(Chem.MolFromSmarts("[C:5](=[CX3H2:6])[C:7](=[O:8])")):
        smarts = "([NX3H1:1][C:2]=[C:3][SX2H1:4].[C:5](=[CX3H2:6])[C:7](=[O:8])) >> [NX3H1:1][C:2]=[C:3][SX2H0:4][CX4H2:6][C:5][C:7](=[O:8])"
        cysquery = 2
    #For labionine formation, needs improvement
    elif bgc_dict["Core Subclass"] == "Class III" and mode == 0:
        smarts = "([C:1]([C:2][SX2H:3])[C:4](=[O:5]).[C:6](=[CX3H2:7])[C:8](=[O:9]).[C:10](=[CX3H2:11])[C:12](=[O:13])) >> [C:1]([C:2][SX2H0:3][CX4H2:7][C:6]([CX4H2:11][C:10][C:12](=[O:13]))[C:8](=[O:9]))[C:4](=[O:5])"
    elif bgc_dict["Core Class"] == "thiopeptide":
        smarts = "([C:1]([C:2][SX2H:3])[C:4](=[O:5]).[C:6](=[CX3H2:7])[N:8][C:9](=[O:10])) >> [C:1]([C:2][SX2H0:3][C:7][C:6][N:8][C:9](=[O:10]))[C:4](=[O:5])"
        dhxquery = 3
        ringsize = 45
    else:
        smarts = "([C:1]([C:2][SX2H:3])[C:4](=[O:5]).[C:6](=[C:7])[C:8](=[O:9])) >> [C:1]([C:2][SX2H0:3][C:7][C:6][C:8](=[O:9]))[C:4](=[O:5])"
    #Creates lists of reacting residues with average position using the remapped values
    mol_serthr = mol.GetSubstructMatches(Chem.MolFromSmarts(queries[dhxquery]))
    mol_cys = mol.GetSubstructMatches(Chem.MolFromSmarts(queries[cysquery]))
    mol_cyslist = average_position(remap(mol_cys, mol_map))
    mol_serthrlist = average_position(remap(mol_serthr, mol_map))

    #Chooses residues best fitting the heuristic (cyclization in C-N direction, from cysteine to ser/thr. 1-2 residues minimum spacing for cycles. NB! 20 (18) is a temporary value, needs to check if this is reasonable)
    #The general way this works is that it picks the cysteine furthest to the C-terminus of the peptide (by using the remapped values, which gives a single average positional value in correct N-C direction), then finds
    #the optimal Dhx by choosing one of the values in the mol_serthrlist which fullfils the requirements, which is that it's closest to the cysteine position while still to the left of it in N-C (plus the value for ringsize)
    if len(mol_cyslist) > 0 and len(mol_serthrlist) > 0 and mol_cyslist[0] < mol_serthrlist[0] and bgc_dict["Core Subclass"] == "Class II":
        optimal = mol_serthrlist[len(mol_serthrlist) - 1]
    else:
        optimal = mol_serthrlist[0]
    rightmost = mol_cyslist[len(mol_cyslist) - 1]
    for st in mol_serthrlist:
        if abs(rightmost - st) < abs(rightmost - optimal) and st < rightmost - ringsize:
            optimal = st

    #Runs reaction
    cyc_rxn = AllChem.ReactionFromSmarts(smarts)
    products = cyc_rxn.RunReactants([mol, ])
    #Creates list containing remapped positional values of reacting cysteine residues for all of the products
    cycmol_cysmaps = []
    for p in products:
        Chem.SanitizeMol(p[0])
        cycmol_cys = p[0].GetSubstructMatches(Chem.MolFromSmarts(queries[cysquery]))
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
        cycmol_serthr = m.GetSubstructMatches(Chem.MolFromSmarts(queries[dhxquery]))
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
    if cyc_rxn.RunReactants([cyc_mol, ]) == () and cysquery == 1 and mode == 0 and bgc_dict["Core Subclass"] == "Class III" and cyc_mol.HasSubstructMatch(Chem.MolFromSmarts(queries[0])) and cyc_mol.HasSubstructMatch(Chem.MolFromSmarts(queries[1])):
        return recursive_cyclization(cyc_mol, n, bgc_dict, 1)
    if cyc_rxn.RunReactants([cyc_mol, ]) == () and cysquery == 1:
        return cyc_mol, n
    if mode == 1:
        return recursive_cyclization(cyc_mol, n, bgc_dict, 1)
    if mode == 0:
        return recursive_cyclization(cyc_mol, n, bgc_dict, 0)

#Function that models the thioether crosslinking reaction between two metabolites. Only supports lanthionine/methyllanthionine formation, labionin/avionin not implemented (completely)
def cyclization(node, model, bgc_dict):
    substrate = node.molecule
    #Counts cycles (pun intended), each corresponds to creation of a water molecule
    n = 0
    product, c = recursive_cyclization(substrate, n, bgc_dict, 0)
    Chem.SanitizeMol(product)
    cyc_smiles = Chem.MolToSmiles(product)

    cyclic_peptide = Metabolite(
        id=node.name+"_cyclic_peptide_c",
        formula=AllChem.CalcMolFormula(product),
        name= node.name+" cyclic peptide",
        compartment="c"
    )

    cyclization_reaction = Reaction("cyclization_of_" + node.name)
    cyclization_reaction.name = "Cyclization reaction of " + node.name
    cyclization_reaction.lower_bound = 0.
    cyclization_reaction.upper_bound = 1000.
    cyclization_reaction.add_metabolites(
        {node.metabolite:-1.0, cyclic_peptide:1.0}
    )

    m_node = Metabolite_node(cyclic_peptide, node.name, node, cyclization_reaction, product, cyc_smiles)

    return m_node

#Models LanO enzymatic tailoring
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
        formula=AllChem.CalcMolFormula(substrate),
        name=node.name + " LanO tailored peptide",
        compartment="c"
    )

    lan_o_tailoring_reaction = Reaction("lan_o_tailoring_reaction_of_"+node.name)
    lan_o_tailoring_reaction.name = "LanO tailoring reaction of "+node.name
    lan_o_tailoring_reaction.lower_bound = 0.
    lan_o_tailoring_reaction.upper_bound = 1000.
    lan_o_tailoring_reaction.add_metabolites(
        {node.metabolite:-1.0, lan_o_tailored_peptide:1.0, UNIMOD.metabolites.get_by_id("h2o_c"):-1.0, UNIMOD.metabolites.get_by_id("nadph_c"):-1.0,
         UNIMOD.metabolites.get_by_id("nadp_c"):1.0, UNIMOD.metabolites.get_by_id("nh4_c"):1.0}
    )

    m_node = Metabolite_node(lan_o_tailored_peptide, node.name, node, lan_o_tailoring_reaction, substrate, lan_o_smiles)
    return m_node

#Models LanD enzymatic tailoring
def lan_d_tailoring(node, model):
    substrate = node.molecule
    cys_thioenol_reaction_smarts = "[NX3H1:1][C:2]([C:3][S:4])[C:5](=[O:6])[OX2H1:7] >> [NX3H1:1][C:2]=[C:3][S:4].[C:5]([O:6])[O:7]"
    ct_reaction = AllChem.ReactionFromSmarts(cys_thioenol_reaction_smarts)
    product = ct_reaction.RunReactants([substrate,])

    Chem.SanitizeMol(product[0][0])
    substrate = product[0][0]

    lan_d_smiles = Chem.MolToSmiles(substrate)

    lan_d_tailored_peptide = Metabolite(
        id=node.name + "_lan_d_tailored_peptide_c",
        formula=AllChem.CalcMolFormula(substrate),
        name=node.name + " LanD tailored peptide",
        compartment="c"
    )

    lan_d_tailoring_reaction = Reaction("lan_d_tailoring_reaction_of_" + node.name)
    lan_d_tailoring_reaction.name = "LanD tailoring reaction of " + node.name
    lan_d_tailoring_reaction.lower_bound = 0.
    lan_d_tailoring_reaction.upper_bound = 1000.
    lan_d_tailoring_reaction.add_metabolites(
        {node.metabolite: -1.0, lan_d_tailored_peptide: 1.0, UNIMOD.metabolites.get_by_id("fmn_c"): -1.0, UNIMOD.metabolites.get_by_id("fmnh2_c"): 1.0,
         UNIMOD.metabolites.get_by_id("co2_c"): 1.0}
    )

    m_node = Metabolite_node(lan_d_tailored_peptide, node.name, node, lan_d_tailoring_reaction, substrate, lan_d_smiles)
    return m_node

#Models thiazole forming heterocyclization reaction
def heterocyclization(node, model, type):
    substrate = node.molecule
    #Counts how many heterocyclizations
    n = 0
    if type == "Type III":
        cysteine_hc_smarts = "[C:1](=[O:2])[N:3][C:4]([C:5][SX2H:6])[C:7][N,n:8] >> [C:1]1=[N:3][C:4](=[C:5][SX2:6]1)[C:7][N,n:8].[O:2]"
    else:
        cysteine_hc_smarts = "[C:1](=[O:2])[N:3][C:4]([C:5][SX2H:6]) >> [C:1]1=[N:3][C:4](=[C:5][SX2:6]1).[O:2]"
    cys_hc_rxn = AllChem.ReactionFromSmarts(cysteine_hc_smarts)
    products = cys_hc_rxn.RunReactants([substrate, ])
    n += 1
    while len(products) > 1:
        Chem.SanitizeMol(products[0][0])
        products = cys_hc_rxn.RunReactants([products[0][0], ])
        Chem.SanitizeMol(products[0][0])
        n += 1

    substrate = products[0][0]

    hc_smiles = Chem.MolToSmiles(substrate)

    thioazole_peptide = Metabolite(
        id=node.name + "_thioazole_containing_peptide_c",
        formula=AllChem.CalcMolFormula(substrate),
        name=node.name + " thioazole containing peptide intermediate",
        compartment="c"
    )
    thioazole_reaction = Reaction("thioazole_heterocyclization_of_" + node.name)
    thioazole_reaction.name = "Thioazole_heterocyclization_reaction_of_" + node.name
    thioazole_reaction.lower_bound = 0.
    thioazole_reaction.upper_bound = 1000.
    thioazole_reaction.add_metabolites(
        {UNIMOD.metabolites.get_by_id("atp_c"): -n, UNIMOD.metabolites.get_by_id("h2o_c"):n, UNIMOD.metabolites.get_by_id("adp_c"): n,
         UNIMOD.metabolites.get_by_id("pi_c"): n, node.metabolite: -1.0, thioazole_peptide: 1.0, UNIMOD.metabolites.get_by_id("fmn_c"): -n,
         UNIMOD.metabolites.get_by_id("fmnh2_c"): n})

    m_node = Metabolite_node(thioazole_peptide, node.name, node, thioazole_reaction, substrate, hc_smiles)
    return m_node

#Reactions that modify the thiopeptide central macrocycle by the Series/Type
def macrocycle_modification(node, model, bgc_dict):
    substrate = node.molecule
    type = bgc_dict["Core Subclass"]
    if type == "Type I":
        smarts = "[N:1][C:2]1[CX4H2:3][C:10]=[C:9]([c:11]3[s:12][c:13][c:14][n:15]3)[N:16]=[C:17]1([c:19]4[n:20][c:21][s:22][c:23]4) >> [C:2]1[CX3H1:3]=[C:10]([O])[C:9]([c:11]3[s:12][c:13][c:14][n:15]3)=[N:16][C:17]=1([c:19]4[n:20][c:21][s:22][c:23]4).[N:1]"
        rxn = AllChem.ReactionFromSmarts(smarts)
        products = rxn.RunReactants([substrate, ])
        product = products[0][0]
        Chem.SanitizeMol(product)
        smiles = Chem.MolToSmiles(product)
        series_e_peptide = Metabolite(
            id=node.name + "_series_e_thiopeptide_c",
            formula=AllChem.CalcMolFormula(product),
            name=node.name + " series e thiopeptide",
            compartment="c"
        )

        modification_reaction = Reaction("Series_E_modification_of_" + node.name)
        modification_reaction.name = "Series E modification reaction of " + node.name
        modification_reaction.lower_bound = 0.
        modification_reaction.upper_bound = 1000.
        modification_reaction.add_metabolites(
            {node.metabolite: -1.0, series_e_peptide: 1.0, UNIMOD.metabolites.get_by_id("h2o_c"): -1.0}
        )
        m_node = Metabolite_node(series_e_peptide, node.name, node, modification_reaction, product, smiles)
        return m_node
    elif type == "Type II":
        smarts = "[N:1][C:2]1[CX4H2:3][C:10]=[C:9]([c:11]3[s:12][c:13][c:14][n:15]3)[N:16]=[C:17]1([c:19]4[n:20][c:21][s:22][c:23]4) >> [N:1][C:2]1[CX4H2:3][C:10][C:9]([c:11]3[s:12][c:13][c:14][n:15]3)=[N:16][C:17]1([c:19]4[n:20][c:21][s:22][c:23]4)"
        rxn = AllChem.ReactionFromSmarts(smarts)
        products = rxn.RunReactants([substrate, ])
        product = products[0][0]
        Chem.SanitizeMol(product)
        smiles = Chem.MolToSmiles(product)
        series_b_peptide = Metabolite(
            id=node.name + "_series_b_thiopeptide_c",
            formula=AllChem.CalcMolFormula(product),
            name=node.name + " series b thiopeptide",
            compartment="c"
        )

        modification_reaction = Reaction("Series_B_modification_of_" + node.name)
        modification_reaction.name = "Series B modification reaction of " + node.name
        modification_reaction.lower_bound = 0.
        modification_reaction.upper_bound = 1000.
        modification_reaction.add_metabolites(
            {node.metabolite: -1.0, series_b_peptide: 1.0, UNIMOD.metabolites.get_by_id("h_c"): -1.0}
        )
        m_node = Metabolite_node(series_b_peptide, node.name, node, modification_reaction, product, smiles)
        return m_node
    elif type == "Type III":
        smarts = "[N:1][C:2]1[CX4H2:3][C:10]=[C:9]([c:11]3[s:12][c:13][c:14][n:15]3)[N:16]=[C:17]1([c:19]4[n:20][c:21][s:22][c:23]4) >> [C:2]1[CX3H1:3]=[C:10][C:9]([c:11]3[s:12][c:13][c:14][n:15]3)=[N:16][C:17]=1([c:19]4[n:20][c:21][s:22][c:23]4).[N:1]"
        rxn = AllChem.ReactionFromSmarts(smarts)
        products = rxn.RunReactants([substrate,])
        product = products[0][0]
        Chem.SanitizeMol(product)
        smiles = Chem.MolToSmiles(product)
        series_d_peptide = Metabolite(
            id=node.name + "_series_d_thiopeptide_c",
            formula=AllChem.CalcMolFormula(product),
            name=node.name + " series d thiopeptide",
            compartment="c"
        )

        modification_reaction = Reaction("Series_D_modification_of_" + node.name)
        modification_reaction.name = "Series D modification reaction of " + node.name
        modification_reaction.lower_bound = 0.
        modification_reaction.upper_bound = 1000.
        modification_reaction.add_metabolites(
            {node.metabolite: -1.0, series_d_peptide: 1.0, UNIMOD.metabolites.get_by_id("h_c"): 1.0}
        )
        m_node = Metabolite_node(series_d_peptide, node.name, node, modification_reaction, product, smiles)
        return m_node
    else:
        return node

#Models macrocyclization reaction of thiopeptides
def macrocyclization(node, model, bgc_dict):    
    substrate = node.molecule
    #SMARTS searching queries, first is for 2 pi N-terminal Dha and second is for 4 pi Dha + thioazole motif
    queries = ["[N:1][C:2](=[CX3H2:3])", "[C:9](=[CX3H2:10])([c:11]2[s:12][c:13][c:14][n:15]2)[N:16][C:17](=[O:18])[c:19]3[n:20][c:21][s:22][c:23]3"]
    # Creates dictionary mapping current molecule idx to old value, with somewhat of the correct N-C direction conserved
    mol_map = map_reacted_residues(substrate, queries)
    sub_4_pi = substrate.GetSubstructMatches(Chem.MolFromSmarts(queries[1]))
    sub_4_pi_list = average_position(remap(sub_4_pi, mol_map))
    sub_2_pi = substrate.GetSubstructMatches(Chem.MolFromSmarts(queries[0]))
    sub_2_pi_list = average_position(remap(sub_2_pi, mol_map))
    optimal_4_pi = sub_4_pi_list[len(sub_4_pi_list)-1]
    optimal_2_pi = sub_2_pi_list[0]
    #Diels alder reaction between 2 and 4 pi components to form a Bycroft-Gowland intermediate
    diels_alder_reaction_smarts = "([N:1][C:2](=[CX3H2:3]).[C:9](=[CX3H2:10])([c:11]2[s:12][c:13][c:14][n:15]2)[N:16][C:17](=[O:18])[c:19]3[n:20][c:21][s:22][c:23]3) >> [N:1][C:2]1[CX4H2:3][C:10]=[C:9]([c:11]3[s:12][c:13][c:14][n:15]3)[N:16]=[C:17]1([c:19]4[n:20][c:21][s:22][c:23]4).[O:18]"
    macro_cyc_rxn = AllChem.ReactionFromSmarts(diels_alder_reaction_smarts)
    products = macro_cyc_rxn.RunReactants([substrate,])
    cycmol_4_pi_maps = []
    cycmol_2_pi_maps = []
    for p in products:
        Chem.SanitizeMol(p[0])
        cycmol_4_pi = p[0].GetSubstructMatches(Chem.MolFromSmarts(queries[1]))
        cycmol_map = map_reacted_residues(p[0], queries)
        cycmol_4_pi_maps.append(average_position(remap(cycmol_4_pi, cycmol_map)))
        cycmol_2_pi = p[0].GetSubstructMatches(Chem.MolFromSmarts(queries[0]))
        cycmol_2_pi_maps.append(average_position(remap(cycmol_2_pi, cycmol_map)))
    #Finds product with least dissimilar positional values to the optimums
    i = 0
    dissim = 0
    for m in range(len(cycmol_4_pi_maps)):
        temp_dissim = 0
        for st in cycmol_4_pi_maps[m]:
            temp_dissim += abs(st - optimal_4_pi)
        for st in cycmol_2_pi_maps[m]:
            temp_dissim += abs(st - optimal_2_pi)
        if temp_dissim > dissim:
            dissim = temp_dissim
            i = m
    cyc_mol = products[i][0]
    cyc_smiles = Chem.MolToSmiles(cyc_mol)

    cyclic_peptide = Metabolite(
        id=node.name + "_macrocyclic_thiopeptide_c",
        formula=AllChem.CalcMolFormula(cyc_mol),
        name=node.name + " macrocyclic thiopeptide",
        compartment="c"
    )

    cyclization_reaction = Reaction("macrocyclization_of_" + node.name)
    cyclization_reaction.name = "Macrocyclization reaction of " + node.name
    cyclization_reaction.lower_bound = 0.
    cyclization_reaction.upper_bound = 1000.
    cyclization_reaction.add_metabolites(
        {node.metabolite: -1.0, cyclic_peptide: 1.0, UNIMOD.metabolites.get_by_id("h2o_c"): 1}
    )

    m_node1 = Metabolite_node(cyclic_peptide, node.name, node, cyclization_reaction, cyc_mol, cyc_smiles)

    m_node2 = macrocycle_modification(m_node1, model, bgc_dict)

    return m_node2

#Sets up necessary pathway of quinaldic acid synthesis for secondary macrocycle addition
def thiostrepton_sec_macrocyc_qa_pathway(node, model):
    qa_1_smiles = "[C:1]1=[C:2][C:3]=[C:4]2[C:5](=[C:6]1)[C:7](=[C:8]([C:16])[N:9]2)[C:10][C:11]([C:12](=[O:13])[O:14])[N:15]"
    temp_mol = Chem.MolFromSmiles(qa_1_smiles)
    AllChem.SanitizeMol(temp_mol)
    qa_metabolite_1 = Metabolite(id=node.name + "_quinaldic_acid_intermediate_1_c", formula= AllChem.CalcMolFormula(temp_mol), compartment="c")
    qa_2_smiles = "[C:1]1=[C:2][C:3]=[C:4]2[C:5](=[C:6]1)[C:7](=[C:8]([C:16])[N:9]2)[C:10][C:11]([C:12](=[O:13])[O:14])=[O:17]"
    temp_mol = Chem.MolFromSmiles(qa_2_smiles)
    AllChem.SanitizeMol(temp_mol)
    qa_metabolite_2 = Metabolite(id=node.name + "_quinaldic_acid_intermediate_2_c", formula=AllChem.CalcMolFormula(temp_mol), compartment="c")
    qa_3_smiles = "[C:1]1=[C:2][C:3]=[C:4][C:5](=[C:6]12)[C:7]([C:8](=[O:19])[C:16])=[C:10][C:11]([C:12](=[O:13])[O:14])=[N:9]2"
    temp_mol = Chem.MolFromSmiles(qa_3_smiles)
    AllChem.SanitizeMol(temp_mol)
    qa_metabolite_3 = Metabolite(id=node.name + "_quinaldic_acid_intermediate_3_c", formula=AllChem.CalcMolFormula(temp_mol), compartment="c")
    qa_4_smiles = "[C:1]1=[C:2][C:3]=[C:4][C:5](=[C:6]12)[C:7]([C:8]([O:19])[C:16])=[C:10][C:11]([C:12](=[O:13])[O:14])=[N:9]2"
    temp_mol = Chem.MolFromSmiles(qa_4_smiles)
    AllChem.SanitizeMol(temp_mol)
    qa_metabolite_4 = Metabolite(id=node.name + "_quinaldic_acid_intermediate_4_c", formula=AllChem.CalcMolFormula(temp_mol), compartment="c")

    qa_reaction_1 = Reaction("quinaldic_acid_intermediate_reaction_1")
    qa_reaction_1.name = node.name+" Quinaldic Acid Intermediate Reaction 1"
    qa_reaction_1.lower_bound = 0.
    qa_reaction_1.upper_bound = 1000.
    qa_reaction_1.add_metabolites(
        {UNIMOD.metabolites.get_by_id("trp__L_c"):-1.0, UNIMOD.metabolites.get_by_id("amet_c"):-1.0, qa_metabolite_1:1.0, UNIMOD.metabolites.get_by_id("ahcys_c"):1.0}
    )
    qa_reaction_2 = Reaction("quinaldic_acid_intermediate_reaction_2")
    qa_reaction_2.name = node.name + " Quinaldic Acid Intermediate Reaction 2"
    qa_reaction_2.lower_bound = 0.
    qa_reaction_2.upper_bound = 1000.
    qa_reaction_2.add_metabolites(
        {qa_metabolite_1:-1.0, qa_metabolite_2:1.0, UNIMOD.metabolites.get_by_id("nh4_c"):1.0}
    )
    qa_reaction_3 = Reaction("quinaldic_acid_intermediate_reaction_3")
    qa_reaction_3.name = node.name + " Quinaldic Acid Intermediate Reaction 3"
    qa_reaction_3.lower_bound = 0.
    qa_reaction_3.upper_bound = 1000.
    qa_reaction_3.add_metabolites(
        {qa_metabolite_2:-1.0, UNIMOD.metabolites.get_by_id("fadh2_c"):-1.0, UNIMOD.metabolites.get_by_id("o2_c"):-1.0,
         qa_metabolite_3:1.0, UNIMOD.metabolites.get_by_id("fad_c"):1.0, UNIMOD.metabolites.get_by_id("h2o_c"):1.0}
    )
    qa_reaction_4 = Reaction("quinaldic_acid_intermediate_reaction_4")
    qa_reaction_4.name = node.name + " Quinaldic Acid Intermediate Reaction 4"
    qa_reaction_4.lower_bound = 0.
    qa_reaction_4.upper_bound = 1000.
    qa_reaction_4.add_metabolites(
        {qa_metabolite_3:-1.0, UNIMOD.metabolites.get_by_id("nadph_c"):-1.0, qa_metabolite_4:1.0, UNIMOD.metabolites.get_by_id("nadp_c"):1.0}
    )
    model.add_reactions([qa_reaction_1, qa_reaction_2, qa_reaction_3, qa_reaction_4])

#Models quinaldic acid addition to thiostrepton-like thiopeptide
def thiostrepton_qa_addition(node, model):
    thiostrepton_sec_macrocyc_qa_pathway(node, model)
    substrate1 = node.molecule
    substrate2 = Chem.MolFromSmarts("[C:1]1=[C:2][C:3]=[C:4][C:5](=[C:6]12)[C:7]([C:8]([O:19])[C:16])=[C:10][C:11]([C:12](=[O:13])[O:14])=[N:9]2")
    Chem.SanitizeMol(substrate1)
    thr_query = ["[C:1]([C:2]([C:3])[OX2H:4])[C:5](=[O:6])"]

    mol_map = map_reacted_residues(substrate1, thr_query)
    #4 binds to 21
    qa_addition_reaction_smarts = "[C:1]([C:2]([C:3])[OX2H:4])[C:5](=[O:6]).[C:7]1=[C:8][C:9]=[C:10][C:11](=[C:12]12)[C:13]([C:14]([O:15])[C:16])=[C:17][C:18]([C:19](=[O:20])[O:21])=[N:22]2 >> [C:7]1=[C:8][C:9]=[C:10][C:11](=[C:12]12)[C:13]([C:14]([O:15])[C:16])=[C:17][C:18]([C:19](=[O:20])[OX2:4]([C:2]([C:3])[C:1][C:5]=[O:6]))=[N:22]2.[O:21]"
    qa_reaction = AllChem.ReactionFromSmarts(qa_addition_reaction_smarts)

    mol_thr = substrate1.GetSubstructMatches(Chem.MolFromSmarts(thr_query[0]))
    mol_thrlist = average_position(remap(mol_thr, mol_map))

    optimum_value = mol_thrlist[len(mol_thrlist)-1]
    products = qa_reaction.RunReactants([substrate1,substrate2])
    product = products[0][0]
    for p in products:
        Chem.SanitizeMol(p[0])
        mol_map = map_reacted_residues(p[0], thr_query)
        mol_thr = p[0].GetSubstructMatches(Chem.MolFromSmarts(thr_query[0]))
        mol_thrlist = average_position(remap(mol_thr, mol_map))
        if optimum_value in mol_thrlist:
            continue
        product = p[0]
        break

    Chem.SanitizeMol(product)
    qa_tpp_smiles = Chem.MolToSmiles(product)

    quinaldic_acid_peptide = Metabolite(
        id=node.name + "_quinaldic_acid_peptide_c",
        formula=AllChem.CalcMolFormula(product),
        name=node.name + " quinaldic acid containing peptide",
        compartment="c"
    )

    quinaldic_acid_reaction = Reaction("quinaldic_acid_addition_reaction_of_"+node.name)
    quinaldic_acid_reaction.lower_bound = 0.
    quinaldic_acid_reaction.upper_bound = 1000.
    quinaldic_acid_reaction.add_metabolites(
        {quinaldic_acid_peptide:1.0, UNIMOD.metabolites.get_by_id("adp_c"):1.0, UNIMOD.metabolites.get_by_id("pi_c"):1.0,
         model.metabolites.get_by_id(node.name + "_quinaldic_acid_intermediate_4_c"):-1.0, UNIMOD.metabolites.get_by_id("atp_c"):-1.0}
    )

    m_node = Metabolite_node(quinaldic_acid_peptide, node.name, node, quinaldic_acid_reaction, product, qa_tpp_smiles)
    return m_node

#Models thiostrepton-like secondary macrocyclization
def thiostrepton_sec_macrocyclization(node, model):
    substrate = node.molecule
    Chem.SanitizeMol(substrate)

    sec_macro_reaction_smarts = "([c:1]1[c:2][c:3][c:4][c:5]([c:6]12)[c:7]([C:8]([O:9])[C:10])[c:11][c:12]([C:13](=[O:14])[O:15])[n:16]2.[NX3H2:17]) >> [C:1]([O])1[C:2]([NX3H1:17])[C:3]=[C:4][C:5](=[C:6]12)[C:7]([C:8]([O:9])[C:10])=[C:11][C:12]([C:13](=[O:14])[O:15])=[N:16]2"
    sec_macro_reaction = AllChem.ReactionFromSmarts(sec_macro_reaction_smarts)

    products = sec_macro_reaction.RunReactants([substrate, ])

    product = products[0][0]
    Chem.SanitizeMol(product)
    sm_pep_smiles = Chem.MolToSmiles(product)

    secondary_macrocycle_peptide = Metabolite(
        id=node.name + "_secondary_macrocycle_peptide_c",
        formula=AllChem.CalcMolFormula(product),
        name=node.name + " secondary macrocycle containing peptide",
        compartment="c"
    )

    secondary_macrocycle_reaction = Reaction("secondary_macrocycle_reaction_of_" + node.name)
    secondary_macrocycle_reaction.lower_bound = 0.
    secondary_macrocycle_reaction.upper_bound = 1000.
    secondary_macrocycle_reaction.add_metabolites(
        {secondary_macrocycle_peptide: 1.0, UNIMOD.metabolites.get_by_id("nadp_c"):1.0, node.metabolite:-1.0, UNIMOD.metabolites.get_by_id("nadph_c"):-1.0,
         UNIMOD.metabolites.get_by_id("o2_c"):-1.0}
    )

    m_node = Metabolite_node(secondary_macrocycle_peptide, node.name, node, secondary_macrocycle_reaction, product, sm_pep_smiles)
    return m_node

#Sets up necessary pathway of 3-methylindolic acid for nosiheptide ring
def nosiheptide_sec_macrocyc_mia_pathway(node, model):
    mia_1_smiles = "[C:1]1=[C:2][C:3]=[C:4]2[C:5](=[C:6]1)[C:7](=[C:8]([C:12](=[O:13])[O:14])[N:9]2)[C:10]"
    temp_mol = Chem.MolFromSmiles(mia_1_smiles)
    AllChem.SanitizeMol(temp_mol)
    mia_metabolite_1 = Metabolite(id=node.name + "_methylindolic_acid_intermediate_1_c",
                                 formula=AllChem.CalcMolFormula(temp_mol), compartment="c")
    mia_2_smiles = "[C:1]1=[C:2][C:3]=[C:4]2[C:5](=[C:6]1)[C:7](=[C:8]([C:12](=[O:13])[S:16])[N:9]2)[C:10]"
    temp_mol = Chem.MolFromSmiles(mia_2_smiles)
    AllChem.SanitizeMol(temp_mol)
    mia_metabolite_2 = Metabolite(id=node.name + "_methylindolic_acid_intermediate_2_c",
                                 formula=AllChem.CalcMolFormula(temp_mol), compartment="c")
    mia_reaction_1 = Reaction("methylindolic_acid_intermediate_reaction_1")
    mia_reaction_1.name = node.name + " Methylindolic Acid Intermediate Reaction 1"
    mia_reaction_1.lower_bound = 0.
    mia_reaction_1.upper_bound = 1000.
    mia_reaction_1.add_metabolites(
        {UNIMOD.metabolites.get_by_id("trp__L_c"): -1.0, mia_metabolite_1: 1.0}
    )
    mia_reaction_2 = Reaction("methylindolic_acid_intermediate_reaction_2")
    mia_reaction_2.name = node.name + " Methylindolic Acid Intermediate Reaction 2"
    mia_reaction_2.lower_bound = 0.
    mia_reaction_2.upper_bound = 1000.
    mia_reaction_2.add_metabolites(
        {mia_metabolite_1: -1.0, mia_metabolite_2: 1.0, UNIMOD.metabolites.get_by_id("adp_c"): 1.0, UNIMOD.metabolites.get_by_id("pi_c"): 1.0,
         UNIMOD.metabolites.get_by_id("h2o_c"): 1.0, UNIMOD.metabolites.get_by_id("atp_c"): -1.0}
    )
    model.add_reactions([mia_reaction_1, mia_reaction_2])

#Models addition of 3-methylindolic acid to peptide
def nosiheptide_mia_addition(node, model, bgc_dict):
    nosiheptide_sec_macrocyc_mia_pathway(node, model)
    substrate1 = node.molecule
    substrate2 = Chem.MolFromSmarts("[C:1]1=[C:2][C:3]=[C:4]2[C:5](=[C:6]1)[C:7](=[C:8]([C:12](=[O:13])[S:16])[N:9]2)[C:10]")
    Chem.SanitizeMol(substrate1)
    cys_query = ["[C:1](=[O:2])[N:3][C:4]([C:5][SX2H:6])"]

    mol_map = map_reacted_residues(substrate1, cys_query)
    mia_addition_reaction_smarts = "[C:1](=[O:2])[N:3][C:4]([C:5][SX2H:6]).[C:7]1=[C:8][C:9]=[C:10]2[C:11](=[C:12]1)[C:13](=[C:14]([C:15](=[O:16])[S:17])[N:18]2)[C:19] >> [C:7]1=[C:8][C:9]=[C:10]2[C:11](=[C:12]1)[C:13](=[C:14]([C:15](=[O:16])[SX2:6][C:5][C:4][N:3][C:1]=[O:2])[N:18]2)[C:19].[S:17]"
    mia_addition_reaction = AllChem.ReactionFromSmarts(mia_addition_reaction_smarts)

    mol_cys = substrate1.GetSubstructMatches(Chem.MolFromSmarts(cys_query[0]))
    mol_cyslist = average_position(remap(mol_cys, mol_map))

    if len(mol_cyslist) > 4:
        optimum_value = mol_cyslist[3]
    else:
        optimum_value = mol_cyslist[len(mol_cyslist)-1]

    products = mia_addition_reaction.RunReactants([substrate1,substrate2])
    product = products[0][0]
    for p in products:
        Chem.SanitizeMol(p[0])
        mol_map = map_reacted_residues(p[0], cys_query)
        mol_cys = p[0].GetSubstructMatches(Chem.MolFromSmarts(cys_query[0]))
        mol_cyslist = average_position(remap(mol_cys, mol_map))
        if optimum_value in mol_cyslist:
            continue
        product = p[0]
        break

    Chem.SanitizeMol(product)
    mia_tpp_smiles = Chem.MolToSmiles(product)

    methylindolic_acid_peptide = Metabolite(
        id=node.name + "_methylindolic_acid_peptide_c",
        formula=AllChem.CalcMolFormula(product),
        name=node.name + " methylindolic acid containing peptide",
        compartment="c"
    )

    methylindolic_acid_reaction = Reaction("methylindolic_acid_addition_reaction_of_" + node.name)
    methylindolic_acid_reaction.lower_bound = 0.
    methylindolic_acid_reaction.upper_bound = 1000.
    methylindolic_acid_reaction.add_metabolites(
        {methylindolic_acid_peptide: 1.0, UNIMOD.metabolites.get_by_id("adp_c"): 1.0,
         UNIMOD.metabolites.get_by_id("pi_c"): 1.0,
         model.metabolites.get_by_id(node.name + "_methylindolic_acid_intermediate_2_c"): -1.0,
         UNIMOD.metabolites.get_by_id("atp_c"): -1.0}
    )

    m_node = Metabolite_node(methylindolic_acid_peptide, node.name, node, methylindolic_acid_reaction, product, mia_tpp_smiles)
    return m_node

#Models nosiheptide-like secondary macrocyclization
def nosiheptide_sec_macrocyclization(node, model):
    #"[c:1]1[c:2][c:3][c:4]2[c:5]([c:6]1)[c:7]([c:8]([C:12](=[O:13])[S:16])[n:9]2)[C:10]"
    substrate = node.molecule
    Chem.SanitizeMol(substrate)
    mia_residue = "[c:7]1[c:8][c:9][c:10]2[c:11]([c:12]1)[c:13]([c:14]([C:15](=[O:16])[SX2:6][C:5][C:4][N:3][C:1]=[O:2])[n:18]2)[C:19]"
    temp_mol = Chem.MolFromSmarts(mia_residue)
    #print(substrate.HasSubstructMatch(temp_mol))

    sec_macro_reaction_smarts = "([CX4H2:1][CX3:2](=[O:3])[OX2H:4].[c:5]1[c:6][c:7][c:8]2[c:9]([c:10]1)[c:11]([c:12]([C:13](=[O:14])[SX2:15][C:16][C:17][N:18][C:19]=[O:20])[n:21]2)[C:22]) >> [c:5]([C][OX2:4][CX3:2](=[O:3])[CX4H2:1])1[c:6][c:7][c:8]2[c:9]([c:10]1)[c:11]([c:12]([C:13](=[O:14])[SX2:15][C:16][C:17][N:18][C:19]=[O:20])[n:21]2)[C:22]"
    sec_macro_reaction = AllChem.ReactionFromSmarts(sec_macro_reaction_smarts)

    products = sec_macro_reaction.RunReactants([substrate,])
    #print(len(products))

    product = products[0][0]
    Chem.SanitizeMol(product)
    sm_pep_smiles = Chem.MolToSmiles(product)

    secondary_macrocycle_peptide = Metabolite(
        id=node.name + "_secondary_macrocycle_peptide_c",
        formula=AllChem.CalcMolFormula(product),
        name=node.name + " secondary macrocycle containing peptide",
        compartment="c"
    )

    secondary_macrocycle_reaction = Reaction("secondary_macrocycle_reaction_of_" + node.name)
    secondary_macrocycle_reaction.lower_bound = 0.
    secondary_macrocycle_reaction.upper_bound = 1000.
    secondary_macrocycle_reaction.add_metabolites(
        {secondary_macrocycle_peptide: 1.0, node.metabolite: -1.0, UNIMOD.metabolites.get_by_id("amet_c"):-1.0, UNIMOD.metabolites.get_by_id("ahcys_c"):1.0}
    )

    m_node = Metabolite_node(secondary_macrocycle_peptide, node.name, node, secondary_macrocycle_reaction, product,
                             sm_pep_smiles)
    return m_node

#Function that finds largest ring substructure in an RDKIT molecule object
def GetLargestRing(mol):
    ssr = Chem.GetSSSR(mol)
    max_ring = 0
    for ring in list(ssr):
        if len(ring) > max_ring:
            max_ring = len(ring)
    return max_ring

#Models macrolactam ring formation for lasso peptides
def macrolactam_ring_formation(node, model):
    substrate = node.molecule

    ml_rxn_smarts = "([CX4H2:1][CX3:2](=[O:3])[OX2H:4].[NX3H2:5][CX4:6]) >> [CX4H2:1][CX3:2](=[O:3])[NX3H1:5][CX4:6].[OX2H2:4]"

    ml_reaction = AllChem.ReactionFromSmarts(ml_rxn_smarts)
    products = ml_reaction.RunReactants([substrate, ])

    maxring = 0
    for p in products:
        if GetLargestRing(p[0]) > maxring:
            product = p[0]
            maxring = GetLargestRing(p[0])

    AllChem.SanitizeMol(product)
    ml_smiles = Chem.MolToSmiles(product)

    macrolactam_peptide = Metabolite(
        id=node.name + "_macrolactam_peptide_c",
        formula=AllChem.CalcMolFormula(product),
        name=node.name + " macrolactam containing peptide",
        compartment="c"
    )

    macrolactam_reaction = Reaction("macrolactam_reaction_of_" + node.name)
    macrolactam_reaction.name = "Macrolactam reaction of " + node.name
    macrolactam_reaction.lower_bound = 0.
    macrolactam_reaction.upper_bound = 1000.
    macrolactam_reaction.add_metabolites(
        {node.metabolite:-1.0, macrolactam_peptide:1.0, UNIMOD.metabolites.get_by_id("atp_c"):-1.0, UNIMOD.metabolites.get_by_id("ppi_c"):1.0,
         UNIMOD.metabolites.get_by_id("amp_c"):1.0}
    )

    m_node = Metabolite_node(macrolactam_peptide, node.name, node, macrolactam_reaction, product, ml_smiles)
    return m_node

#Function for disulfide bridging of class I, III and IV lasso peptides
def disulfide_bridging(node, model):
    substrate = node.molecule
    Chem.SanitizeMol(substrate)
    cysquery = ["[C:1]([C:2][SX2H:3])[C:4](=[O:5])"]

    mol_map = map_reacted_residues(substrate, cysquery)
    ds_smarts = "([SX2H:1].[SX2H:2]) >> [SX2:1][SX2:2]"
    ds_reaction = AllChem.ReactionFromSmarts(ds_smarts)

    mol_cys = substrate.GetSubstructMatches(Chem.MolFromSmarts(cysquery[0]))
    mol_cyslist = average_position(remap(mol_cys, mol_map))
    if len(mol_cyslist) > 2:
        n=2.0
        optimum_coordinating_value = mol_cyslist[0]+mol_cyslist[2]
        products = ds_reaction.RunReactants([substrate,])
        for p in products:
            Chem.SanitizeMol(p[0])
            mol_map = map_reacted_residues(p[0], cysquery)
            mol_cys = p[0].GetSubstructMatches(Chem.MolFromSmarts(cysquery[0]))
            mol_cyslist = average_position(remap(mol_cys, mol_map))
            if mol_cyslist[0]+mol_cyslist[1] == optimum_coordinating_value:
                product = p[0]
                break
            else:
                product = p[0]
        products = ds_reaction.RunReactants([product,])
        product = products[0][0]
    else:
        n=1.0
        products = ds_reaction.RunReactants([substrate, ])
        product = products[0][0]

    AllChem.SanitizeMol(product)
    ds_smiles = Chem.MolToSmiles(product)

    disulfide_bridged_peptide = Metabolite(
        id=node.name + "_disulfide_bridged_peptide_c",
        formula=AllChem.CalcMolFormula(product),
        name=node.name + " disulfide bridge containing peptide",
        compartment="c"
    )

    disulfide_bridge_reaction = Reaction("disulfide_bridge_reaction_of_" + node.name)
    disulfide_bridge_reaction.name = "Disulfide bridge reaction of " + node.name
    disulfide_bridge_reaction.lower_bound = 0.
    disulfide_bridge_reaction.upper_bound = 1000.
    disulfide_bridge_reaction.add_metabolites(
        {node.metabolite: -1.0, disulfide_bridged_peptide: 1.0, UNIMOD.metabolites.get_by_id("gthox_c"):-n, UNIMOD.metabolites.get_by_id("gthrd_c"):2*n}
    )

    m_node = Metabolite_node(disulfide_bridged_peptide, node.name, node, disulfide_bridge_reaction, product, ds_smiles)
    return m_node

#Reaction to indicate the mature product
def maturation(node, model, bgc_dict):
    substrate = node.molecule
    AllChem.SanitizeMol(substrate)
    m_smiles = Chem.MolToSmiles(substrate)
    mature_peptide = Metabolite(
        id=bgc_dict["Name"]+"_"+bgc_dict["Core Number"]+"_mature_peptide_c",
        name="mature peptide",
        formula= AllChem.CalcMolFormula(substrate),
        compartment="c"
    )

    mature_reaction = Reaction(bgc_dict["Name"]+"_"+bgc_dict["Core Number"]+"_maturation_reaction")
    mature_reaction.lower_bound = 0.0
    mature_reaction.upper_bound = 1000.0
    mature_reaction.add_metabolites(
        {node.metabolite: -1.0, mature_peptide: 1.0})

    m_node = Metabolite_node(mature_peptide, node.name, node, mature_reaction, substrate, m_smiles)
    return m_node

#Returns tanimoto similarity for the fingerprints of two SMILES structures
def match(predicted_smiles, BGC_ID):
    mols = [Chem.MolFromSmiles(predicted_smiles), Chem.MolFromSmiles(BGC_to_smiles[BGC_ID])]
    fps = [Chem.RDKFingerprint(x) for x in mols]
    return DataStructs.FingerprintSimilarity(fps[0], fps[1])

#Returns list of records from .gbk file (taken from original BiGMeC script)
def get_gb_list_from_antismash_output(cluster_path):  # yes
    # domains, pfam_entries, smCOG, EC_number, rxn_in, rxn_out, core_gene, amino_acid_sequence, predicted_EC, length
    # triple pound signs means that information is stored within a custom class object.
    gb_list = []
    for gb_record in SeqIO.parse(open(cluster_path, "r"), "genbank"):
        gb_list.append(gb_record)
    return gb_list

#Checks for presence of RiPP prepeptide cores in .gbk record features, runs pathway reconstruction for each
def analyse_bgc(bgc_path, bgc_dict, output, mode):
    result_list = []
    gb_list, bgc_dict = parse_gbk(bgc_path, bgc_dict)
    cores = 1
    #Runs for each core identified in the BGC file, with same conditions
    for gb_record in gb_list:
        for feat in gb_record.features:
            if feat.type == "CDS_motif":
                if feat.qualifiers.get("prepeptide") == ["core"] and cores == 1:
                    bgc_dict["Core"] = feat.qualifiers.get("core_sequence")
                    bgc_dict["Leader"] = feat.qualifiers.get("leader_sequence")
                    bgc_dict["Core Number"] = str(cores)
                    bgc_dict["Core Name"] = feat.qualifiers.get("locus_tag")[0]
                    bgc_dict["Core Subclass"] = feat.qualifiers.get("predicted_class")[0]
                    bgc_dict["Core Class"] = feat.qualifiers.get("peptide")[0]
                    cores +=1
                    result = dict(_run(output, bgc_dict, mode))
                    result_list.append(result)
                elif feat.qualifiers.get("prepeptide") == ["core"]:
                    bgc_dict["Core"] = feat.qualifiers.get("core_sequence")
                    bgc_dict["Leader"] = feat.qualifiers.get("leader_sequence")
                    bgc_dict["Core Number"] = str(cores)
                    bgc_dict["Core Name"] = feat.qualifiers.get("locus_tag")[0]
                    bgc_dict["Core Subclass"] = feat.qualifiers.get("predicted_class")[0]
                    bgc_dict["Core Class"] = feat.qualifiers.get("peptide")[0]
                    cores += 1
                    result = dict(_run(output, bgc_dict, mode))
                    result_list.append(result)
    return result_list

#Depending on run mode, either chooses appropriate CarveMe model associated with BGC source organism if present, otherwise uses Sco-GEM reference model
def model_select(mode, bgc_dict):
    if mode == 1:
        model = cobra.io.read_sbml_model("Models\\Sco-GEM.xml")
        model_type = "REF"
        return model, model_type
    path = "Models\\CarveMe"
    model = None
    model_type = "Unknown"
    tax_id = "Unknown"
    acc_id = "Unknown"
    if bgc_dict["Name"] in TAXID_DICT.keys():
        tax_id = TAXID_DICT[bgc_dict["Name"]]
    if tax_id in NCBI_ACCID_DICT.keys():
        acc_id = NCBI_ACCID_DICT[tax_id]
    for filename in os.listdir(path):
        if acc_id in filename:
            model = cobra.io.read_sbml_model(os.path.join(path, filename))
            model_type = "SRC"
            break
    if model == None:
        model = cobra.io.read_sbml_model("Models\\Sco-GEM.xml")
        model_type = "REF"
    return model, model_type

#Main run code, sets up BGC dictionary and summary dataframe. Runs reconstruction for every .gbk file if ran for folder.
def run(bgc_path, output, mode):
    print("Running...")
    summary = pd.DataFrame({"Name":[], "Core Name":[], "Core Class":[], "Core Subclass":[], "Core Number":[], "BGC":[], "Organism":[], "End product SMILES":[], "Metabolite list":[], "Prediction Score":[], "Precursor Length":[]})
    column_list = ["Name", "Core Name", "Core Class", "Core Subclass", "Core Number", "BGC", "Organism", "End product SMILES", "Metabolite list", "Base Prediction", "Prediction Score", "Precursor Length"]
    bgc_path = Path(bgc_path)
    result_list = []

    if bgc_path.is_dir():
        for filename in tqdm(os.listdir(bgc_path)):
            p = os.path.join(bgc_path, filename)
            if os.path.isfile(p):
                #print(filename)
                bgc_dict = {
                "Name":os.path.splitext(filename)[0], "Core Name":None, "Core Number":None, "BGC":"Unknown", "Organism":"Unknown", "End product SMILES":"Unknown", "Metabolite list":[],
                "SMCOG1070":False, "SMCOG1140":False, "SMCOG1155":False, "SMCOG1075":False, "SMCOG1030":False, "Core":None, "Leader":None, "LanO":False, "Base Prediction":None, "Prediction Score":None,
                "LanD":False, "Core Class":"Unknown", "Core Subclass":"Unknown", "YcaO":False, "Oxazole":False, "Thiostrepton":False, "Precursor Length":0
                }
                result_list.extend(analyse_bgc(p, bgc_dict, output, mode))
    else:
        name = os.path.basename(bgc_path)
        bgc_dict = {
            "Name": os.path.splitext(name)[0], "Core Number":None, "Core Name":None, "BGC": "Unknown", "Organism": "Unknown", "Init product SMILES":"Unknown",
            "End product SMILES": "Unknown", "Metabolite list": [], "SMCOG1070": False, "SMCOG1140": False, "SMCOG1155": False, "SMCOG1075": False, "SMCOG1030":False,
            "Core": None, "Leader": None, "LanO":False, "Base Prediction":None, "Prediction Score":None,
            "LanD":False, "Core Class":"Unknown", "Core Subclass":"Unknown", "YcaO":False, "Oxazole":False, "Thiostrepton":False
        }
        result_list.extend(analyse_bgc(bgc_path, bgc_dict, output, mode))

    for result in result_list:
        new_row = {}
        for column in column_list:
            new_row[column] = [result[column]]
        temp_df = pd.DataFrame(new_row)
        summary = pd.concat([temp_df, summary.loc[:]]).reset_index(drop=True)

    summary.to_csv(os.path.join(output, "RiPP reconstruction summary.csv"))

#Iterates through .gbk file features and updates BGC dictionary (reaction environment). Partially taken from the main BiGMec script
def parse_gbk(bgc_path, bgc_dict):
    count = 0
    gb_list = get_gb_list_from_antismash_output(bgc_path)
    for gb_record in gb_list:
        bgc_dict["Organism"] = gb_record.annotations["organism"]
        for feat in gb_record.features:
            if feat.type == "CDS":
                if feat.qualifiers.get("gene_functions") != None:
                    for i in feat.qualifiers.get("gene_functions"):
                        if "thiopeptide: thiostrepton" in i:
                            bgc_dict["Thiostrepton"] = True
                        if "thiopeptide: YcaO" in i:
                            bgc_dict["YcaO"] = True
                        if "thiopeptide: TIGR03604" in i:
                            count += 1
                        if "smcogs" in i:
                            pos = i.find("smcogs") + 8
                            bgc_dict[i[pos:pos + 9]] = True
                            if feat.qualifiers.get("gene") != None and i[pos:pos + 9] in ["SMCOG1007", "SMCOG1001", "SMCOG1251"] and feat.qualifiers.get("gene")[0][3:4] == "O":
                                bgc_dict["LanO"] = True
                if feat.qualifiers.get("gene") != None and feat.qualifiers.get("gene")[0][3:4] == "D":
                    bgc_dict["LanD"] = True
            elif feat.type == "misc_feature":
                bgc_dict["BGC"] = feat.qualifiers.get("note")
    if count == 2:
        bgc_dict["Oxazole"] = True
    return gb_list, bgc_dict

#Run code for every identified core in the .gbk file. Runs for every CarveMe model if mode 3 is selected
def _run(output, bgc_dict, mode):
    if bgc_dict["Core"] == None and bgc_dict["Leader"] == None:
        bgc_dict["Core Name"] = "Missing core and leader peptide sequence"
        return bgc_dict
    elif bgc_dict["Core"] == None:
        bgc_dict["Core Name"] = "Missing core peptide sequence"
        return bgc_dict
    elif bgc_dict["Leader"] == None:
        bgc_dict["Core Name"] = "Missing leader peptide sequence"
        return bgc_dict
    #Run main script here using altered dictionary, main script saves sbml and PIL files in same output folder
    if mode == 3:
        path = "Models\\CarveMe"
        for m in os.listdir(path):
            #print(m)
            model = cobra.io.read_sbml_model(os.path.join(path, m))
            try:
                bgc_dict = add_to_model(bgc_dict, output, mode, model, os.path.splitext(m)[0])
            except:
                print("Error")
    else:
        bgc_dict = add_to_model(bgc_dict, output, mode, None, None)
    return bgc_dict

#Main pathway reconstruction part of script. Chooses model based on mode, creates chain of metabolite nodes using the reaction functions.
#Iterates through metabolite node chain adding reactions to the model object, returns pathway extended model to output folder
def add_to_model(bgc_dict, output, mode, hetmodel, modelname):
    if mode == 3:
        model = hetmodel
        model_type = "HET"
    elif mode == 4:
        model = Model(bgc_dict["Name"]+"_"+bgc_dict["Core Number"])
        model_type = "RPW"
    else:
        model, model_type = model_select(mode, bgc_dict)
    leaderseq = bgc_dict["Leader"][0].replace(" ","").replace("X","G")
    coreseq = bgc_dict["Core"][0].replace(" ","").replace("X","G")
    pepseq = leaderseq+coreseq
    bgc_dict["Precursor Length"] = len(pepseq)

    origin_node = translation(pepseq, coreseq, bgc_dict["Name"], model)
    current_node = origin_node
    #Quinaldic acid addition for later secondary macrocyclization of thiostrepton-like thiopeptides
    if bgc_dict["Core Class"] == "thiopeptide" and bgc_dict["Core Subclass"] == "Type II" and bgc_dict["Thiostrepton"]:
        current_node = thiostrepton_qa_addition(current_node, model)
    #Methylindolic acid addition and secondary macrocyclization of nosiheptide-like thiopeptides
    if bgc_dict["Core Class"] == "thiopeptide" and bgc_dict["Core Subclass"] == "Type I" and bgc_dict["Thiostrepton"]:
        current_node = nosiheptide_mia_addition(current_node, model, bgc_dict)
        current_node = nosiheptide_sec_macrocyclization(current_node, model)
    #Heterocyclization of cysteine residues in thiopeptides (producing thioazoles, oxazoles not implemented)
    if heterocyclizationable(current_node) and bgc_dict["YcaO"] or heterocyclizationable(current_node) and bgc_dict["Core Class"] == "thiopeptide":
        current_node = heterocyclization(current_node, model, bgc_dict["Core Subclass"])
    #Glutamylation dehydration of either class I lanthipeptides or thiopeptides
    if dehydratable(current_node) and bgc_dict["SMCOG1155"] or dehydratable(current_node) and bgc_dict["Core Class"] == "thiopeptide":
        current_node = glutamylation_elimination(current_node, model, 1)
    #Phosphorylation dehdyration of class II, III or IV lanthipeptides
    if dehydratable(current_node) and bgc_dict["SMCOG1030"] or dehydratable(current_node) and bgc_dict["SMCOG1070"]:
        current_node = phosphorylation_dehydration(current_node, model, 1)
    #LanD mediated C-terminus tailoring of class I lanthipeptides
    if c_term_oxidative_decarboxylationable(current_node) and bgc_dict["LanD"] and bgc_dict["Core Class"] == "lanthipeptide":
        current_node = lan_d_tailoring(current_node, model)
    #Cyclization of lanthipeptides
    if cyclizationable(current_node) and bgc_dict["SMCOG1140"] or cyclizationable(current_node) and bgc_dict["SMCOG1070"] or cyclizationable(current_node) \
            and bgc_dict["SMCOG1030"] or cyclizationable(current_node) and bgc_dict["Core"] == "thiopeptide" and bgc_dict["Core Subclass"] == "Type III":
        current_node = cyclization(current_node, model, bgc_dict)
    #LanO mediated N-terminus tailoring of class I lanthipeptides
    if n_term_lac_able(current_node) and bgc_dict["LanO"] and bgc_dict["Core Class"] == "lanthipeptide":
        current_node = lan_o_tailoring(current_node, model)
    #Primary macrocyclization of thiopeptide
    if macrocyclizationable(current_node) and bgc_dict["Core Class"] == "thiopeptide":
        current_node = macrocyclization(current_node, model, bgc_dict)
    #Secondary macrocyclization of thiostrepton-like thiopeptides
    if bgc_dict["Core Class"] == "thiopeptide" and bgc_dict["Core Subclass"] == "Type II" and bgc_dict["Thiostrepton"]:
        current_node = thiostrepton_sec_macrocyclization(current_node, model)
    #Macrolactam reaction of lassopeptides
    if macrolactamable(current_node) and bgc_dict["Core Class"] == "lassopeptide":
        current_node = macrolactam_ring_formation(current_node, model)
    #Disulfide bridging of class I, III and IV lasso peptides
    if disulfide_bridgable(current_node) and bgc_dict["Core Class"] == "lassopeptide":
        current_node = disulfide_bridging(current_node, model)
    if bgc_dict["Core Class"] == "lanthipeptide" or bgc_dict["Core Class"] == "lassopeptide":
        current_node = core_cleave(current_node, leaderseq, model, bgc_dict)
    #Need to move core/leader cleavage here, differentiate by class
    current_node = maturation(current_node, model, bgc_dict)
    #current_node = transport(current_node, model)
    bgc_dict["End product SMILES"] = current_node.smiles

    if mode == 1 and bgc_dict["Name"] in BGC_to_smiles.keys():
        bgc_dict["Prediction Score"] = match(bgc_dict["End product SMILES"], bgc_dict["Name"])
        #print(bgc_dict["Prediction Score"])
    if mode == 2:
        m = 0
        for k in BGC_to_smiles.keys():
            if match(bgc_dict["End product SMILES"], k) > m:
                m = match(bgc_dict["End product SMILES"], k)
                bgc_dict["Prediction Score"] = k

    metabolite_list = []
    while current_node != None:
        model.add_reactions([current_node.reaction])
        metabolite_list.append(current_node.metabolite.id)
        if current_node.parent == None and mode == 1 and bgc_dict["Name"] in BGC_to_smiles.keys():
            bgc_dict["Base Prediction"] = match(current_node.smiles, bgc_dict["Name"])
        current_node = current_node.parent
    metabolite_list.reverse()

    bgc_dict["Metabolite list"] = metabolite_list
    bgc_dict["Success"] = "Yes"

    #Creates mature peptide demand reaction for flux analysis
    model.add_boundary(model.metabolites.get_by_id(metabolite_list[len(metabolite_list)-1]), type="demand")
    #Saves model for each core
    if mode == 3:
        cobra.io.write_sbml_model(
            model, join(output, "Pathway Models\\" + model_type + "_" + bgc_dict["Name"] + "_" + bgc_dict["Core Number"] + "_" + modelname + ".xml")
        )
    else:
        cobra.io.write_sbml_model(model, join(output,"Pathway Models\\"+model_type+"_"+bgc_dict["Name"]+"_"+bgc_dict["Core Number"]+".xml"))
    return bgc_dict

def main():
    #Mode 1: Matches both end product and precursor peptide with characterized structures, uses only reference model
    #Mode 2: Scramble test, matches with the best fitting of the set of characterized structures. Returns the associated BGC ID instead of score
    #Mode 3: Heterologous expression extension for all CarveMe models
    #Mode 4: Exports the generated pathways as a set of minimodels only containing the reactions, for later use
    run("Data\\MIBiG RiPP gbks",
        "RiPP Pathway Output", 4)

if __name__ == "__main__":
    main()