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

BGC_to_smiles = {
    "BGC0002111":r"C/C=C(\NC(=O)/C(=C\C)NC(=O)CNC(=O)C(NC(=O)C(CC1=C[NH]C2=CC=CC=C12)NC(=O)C(C)NC1CC2(CSCC(C(=O)O)NC(=O)C(CC(N)=O)NC(=O)C(CC(N)=O)NC2=O)NC(=O)C(C(C)CC)NC(=O)/C(=C/C)NC1=O)C(C)C)C(=O)N/C(=C/C)C(=O)NC(C(=O)N/C(=C/C)C(=O)NC(CC1=C[NH]C2=CC=CC=C12)C(=O)N/C(=C/C)C(=O)N/C(=C\C)C(=O)N/C(=C/C)C(=O)NC(C(=O)N/C(=C/C)C(=O)NC(C)C(=O)N(C)C)C(C)CC)C(C)C",
    "BGC0001579":r"C[C@@H]1[C@H]2C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H]3CSC[C@H]4C(=O)N[C@@H](CS1)C(=O)N[C@@H](CNCCCC[C@H](NC(=O)[C@H]([C@H](SC[C@@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N4)CCC(=O)N)CCCCN)N)C)NC(=O)[C@@H](NC(=O)CNC(=O)[C@@H](NC3=O)[C@H](C(=O)O)O)CC(=O)N)C(=O)O)C(=O)N[C@H](C(=O)NCC(=O)N5CCC[C@H]5C(=O)N[C@H](C(=O)N2)CC6=CC=CC=C6)CC7=CC=CC=C7)C(C)C)CC8=CC=CC=C8",
    "BGC0000559":r"NC(CC1=CNC2=C1C=CC=C2)C(NC(CCCCN)C(NC(CSCC(C(NC(C5C)C(N4C(C(NC([H])C(NC(CS5)C(NC(C(C)C)C(NC(C6C)C(NC([H])C(NC(C)C(NC(CC(C)C)C(NC(CCC(N)=O)C(NC(C(NC(CS6)C(NC(CC7=CC=CC=C7)C(NC(CC(C)C)C(NC(CCC(N)=O)C(NC(C8C)C(NC(CC(C)C)C(NC(C9C)C(NC(CS8)C(NC(CC(N)=O)C(NC(CS9)C(NC(CCCCN)C(NC(C(CC)C)C(NC(C(NC(CCCCN)C(O)=O)=O)=C)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=CC)=O)=O)=O)=O)=O)=O)=O)=O)=O)CCC4)=O)=O)NC3=O)C(NC(CCC(O)=O)C(NC(C(NC3CC(C)C)=O)=C)=O)=O)=O)=O",
    "BGC0000551":r"CCC(C)C1C(=O)NC(C(=O)NC(C(=O)NC(CSCC(C(=O)NC(C(=O)NC(C(=O)NC(=C)C(=O)N1)CC(C)C)CO)NC(=O)C(CC(=O)O)NC(=O)CNC(=O)C2CSCC(C(=O)NC(C(=O)NC(C(=O)NC(=C)C(=O)NC(C(=O)NC(C(=O)NC(C(=O)N2)CC(C)C)CC(C)C)CC(C)C)C)CCCNC(=N)N)NC(=O)CNC(=O)C(C(C)O)N)C(=O)NC(CC(=O)N)C(=O)O)C(C)O)C(C)O",
    "BGC0000535":r"CCC(C)C1C(=O)NC(=C)C(=O)NC(C(=O)NC(CSCC(C(=O)N1)NC(=O)/C(=C/C)/NC(=O)C(C(C)CC)N)C(=O)NC2C(SCC(NC(=O)CNC(=O)C3CCCN3C2=O)C(=O)NC(CCCCN)C(=O)NC4C(SCC(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC4=O)C)CC(C)C)CCSC)C(=O)NC(CC(=O)N)C(=O)NC(CCSC)C(=O)NC(CCCCN)C(=O)NC5C(SCC6C(=O)NC(C(=O)NC(CSC(C(C(=O)N6)NC(=O)C(NC5=O)C)C)C(=O)NC(CO)C(=O)NC(C(C)CC)C(=O)NC(CC7=CNC=N7)C(=O)NC(C(C)C)C(=O)NC(=C)C(=O)NC(CCCCN)C(=O)O)CC8=CNC=N8)C)C)C)CC(C)C",
    "BGC0000528":r"CC[C@H](C)[C@@H]1NC(=O)[C@H](CCSC[C@@H]2NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](CCSC[C@H](NC(=O)[C@H](C)NC2=O)C(=O)N[C@@H](CCCNC(N)=N)C(O)=O)NC(=O)CNC(=O)[C@H](CS)NC(=O)[C@H](CCC(O)=O)NC1=O)[C@@H](C)CC)[C@@H](C)CC)NC(=O)[C@H](CC(C)C)NC(=O)C(\NC(=O)[C@@H]1CSC[C@H](NC(=O)[C@@H](N)CO)C(=O)N[C@@H](CO)C(=O)NCC(=O)N[C@@H](CC2=CNC3=C2C=CC=C3)C(=O)N[C@@H](CC(C)C)C(=O)N1)=C\C",
    "BGC0000523":r"CSCCC(NC(=O)C(C)NC(=O)C(NC(=O)C(C)NC(=O)C(NC(=O)C(C)NC(=O)C(C)NC(=O)C(CC(C)C)NC(=O)C(NC(=O)C1CCCN1C(=O)C(\NC(=O)C(C)=O)=C/C)C(C)C)C(C)C)C(C)C)C(=O)NC(CCC(O)=O)C(=O)NC(CC(C)C)C(=O)NC(CC(C)C)C(=O)N1CCCC1C(=O)NC(C(C)O)C(=O)NC(C)C(=O)NC(C)C(=O)NC(C(C)C)C(=O)NC(CC(C)C)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC1CSCC(NC(=O)CNC(=O)C(C)NC(=O)C(NC(=O)C(CC(O)=O)NC1=O)C(C)C)C(=O)NC(CC1=CC=CC=C1)C(=O)NC(CCCCN)C(=O)NC(CC1=CC=C(O)C=C1)C(=O)NC1CSCC(NC(=O)C(CC2=CNC=N2)NC(=O)C(CC2=CNC=N2)NC(=O)C(CCCCN)NC(=O)C(C)NC1=O)C(O)=O",
    "BGC0000521":r"CCC(C)C(NC(=O)C(NC(=O)CNC(=O)C(CO)NC(=O)CNC(=O)CNC(=O)C(N)CCCCN)C(C)C)C(=O)NC(CC1=CNC=N1)C(=O)NC1C(C)SCC2NC(=O)C(CCC(O)=O)NC(=O)C(CC3=CNC=N3)NC(=O)C(CSCC3NC(=O)\C(NC(=O)C(CC4=CC=CC=C4)NC(=O)C(NC(=O)C(CC4=CC=CC=C4)NC(=O)C(CCC(N)=O)NC(=O)C(CC4=CNC5=C4C=CC=C5)NC(=O)C(CSCC(NC3=O)C(=O)NC(CO)C(O)=O)NC(=O)C(CC(N)=O)NC(=O)C(CCSC)NC(=O)C(CC(N)=O)NC2=O)C(C)C)=C\C)NC(=O)C(NC1=O)C(C)CC",
    "BGC0000511":r"NC(CC1=CNC2=C1C=CC=C2)C(NC(CCCCN)C(NC(CSCC(C(NC(C5C)C(N4C(C(NC([H])C(NC(CS5)C(NC(C(C)C)C(NC(C6C)C(NC([H])C(NC(C(C)C)C(NC(CC(C)C)C(NC(CCC(N)=O)C(NC(C(NC(CS6)C(NC(CC7=CC=CC=C7)C(NC(CC(C)C)C(NC(CCC(N)=O)C(NC(C8C)C(NC(C(CC)C)C(NC(C9C)C(NC(CS8)C(NC(CC(N)=O)C(NC(CS9)C(NC(CC%10=CNC=N%10)C(NC(C(CC)C)C(NC(C(NC(CCCCN)C(O)=O)=O)=C)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=CC)=O)=O)=O)=O)=O)=O)=O)=O)=O)CCC4)=O)=O)NC3=O)C(NC(CCC(O)=O)C(NC(C(NC3C(C)C)=O)=C)=O)=O)=O)=O",
    "BGC0000509":r"CC[C@H](C)[C@H](NC(=O)C(\NC(=O)C(\NC(=O)[C@H](CCCCN)NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)C(=C)NC(=O)[C@H](C)NC(=O)C(C)O)[C@@H](C)CC)C(C)C)=C/C)=C/C)C(=O)NC(CCCCN)C(=O)N[C@@H](C)C(=O)N[C@H]1CSC[C@H](NC(=O)[C@H](CC(C)C)NC(=O)C(CCCCN)NC(=O)C(CCCCN)NC1=O)C(=O)N[C@@H](CCCNC(N)=N)C(=O)NCC(=O)N[C@@H](CC1=CC=CC=C1)C(=O)N[C@H]1C(C)SC[C@@H]2NC(=O)[C@@H](NC(=O)[C@H](CC(C)C)NC1=O)C(C)SC[C@H](NC(=O)CNC2=O)C(=O)N[C@@H](CC1=CN=CN1)C(=O)N[C@@H](CC1=CC=CC=C1)C(=O)N\C(=C\C)C(=O)NCC(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(O)=O",
    "BGC0000508":r"CCC(C)C1C(=O)NC(CSCC(C(=O)NC(C(=O)NC(C(=O)N1)CC2=CC=CC=C2)CCCCN)NC(=O)C(C)NC(=O)C(C(C)CC)N)C(=O)NC3C(SCC(NC(=O)CNC(=O)C4CCCN4C3=O)C(=O)NC(C)C(=O)NC(CCCCN)C(=O)N/C(=C\C)/C(=O)NCC(=O)NC5CSCC6C(=O)N/C=C\SCC(C(=O)NC(C(=O)N6)CC7=CC=C(C=C7)O)NC(=O)C(NC(=O)C(NC5=O)CC8=CC=CC=C8)CC(=O)N)C",
    "BGC0000503":r"C[C@H]1[C@@H]2C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H]3CSC[C@@H]4C(=O)N[C@@H](CS1)C(=O)N[C@@H](CNCCCC[C@H](NC(=O)[C@@H]([C@@H](SC[C@@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N4)CCC(=O)N)CCCNC(=N)N)N)C)NC(=O)[C@@H](NC(=O)CNC(=O)[C@@H](NC3=O)[C@@H](C(=O)O)O)CC(=O)N)C(=O)O)C(=O)N[C@H](C(=O)NCC(=O)N5CCC[C@H]5C(=O)N[C@H](C(=O)N2)CC6=CC=CC=C6)CC7=CC=CC=C7)C(C)C)CC8=CC=CC=C8",
    "BGC0000501":r"CC[C@H](C)[C@@H]1NC(=O)CNC(=O)[C@H](C[C@@]2(CSC[C@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(C)C)NC2=O)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)NC1=O)NC(=O)[C@H](CC1=CN=CN1)NC(=O)CNC(=O)[C@H](CC(N)=O)NC(=O)[C@@H]1CSC[C@@]2(C[C@H](NC(=O)[C@H](CC(O)=O)NC(=O)CNC(=O)CNC(=O)CNC(=O)CNC(=O)[C@H](CC3=CN=CN3)NC(=O)CN)C(=O)NCC(=O)N[C@@H](CC(C)C)C(=O)N2)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H]([C@@H](C)O)C(=O)NCC(=O)N1",
    "BGC0000549":r"[H]N[C@@H](CC1=CC=CC=C1)C(N/C(C(N[C@@H]2C(N[C@@H](CC3=CN=CN3)C(NC(C(N[C@H](C(N[C@@H](CSC2)C(N[C@@H]4C(N5CCC[C@H]5C(NCC(N[C@@H](CSCC4)C(N[C@]([C@@H](C)CC)([H])C(N[C@@H]6C(NCC(N[C@H](C(N[C@H](C(N[C@@H](CCSC)C(NCC(N[C@@H](CSCC6)C(N[C@@H](CC7=CN=CN7)C(N[C@]([C@@H](C)CC)([H])C(N[C@H](C(N[C@@H]8C(N[C@]([C@@H](C)CC)([H])C(NCC(N[C@@H](CSC8)C(N[C@H](C(N[C@H](C(N[C@@H](CC9=CN=CN9)C(N[C@]([C@@H](C)CC)([H])C(N[C@H](C(N[C@]([H])(C(NC(C(N[C@H](C(O)=O)CCCCN)=O)=C)=O)[C@H](CC)C)=O)CC%10=CN=CN%10)=O)=O)=O)C(C)C)=O)CC(N)=O)=O)=O)=O)=O)=O)CCC(N)=O)=O)=O)=O)=O)=O)=O)CC(C)C)=O)C(C)C)=O)=O)=O)=O)=O)=O)=O)=O)=O)CC(C)C)=O)=C)=O)=O)=O)=C\C)=O",
    "BGC0000544":r"CC[C@H](C)[C@H](N)C(=O)N\C(=C/C)C(=O)N[C@H]1CSC[C@@H]2NC(=O)[C@@H](NC(=O)[C@@H]3CSCC[C@H](NC(=O)[C@H](CSC[C@@H]4NC(=O)CNC(=O)[C@H](CSC[C@@H]5NC(=O)[C@H](CC6=CN=CN6)NC(=O)[C@H](CSC[C@H](NC5=O)C(=O)N5CCC[C@H]5C(=O)NCC(O)=O)NC4=O)NC(=O)CNC(=O)CNC(=O)CNC(=O)[C@H](CCC(O)=O)NC2=O)NC(=O)[C@H](CC2=CNC4=C2C=CC=C4)NC(=O)C(=C)NC(=O)[C@@H](NC1=O)C(C)C)C(=O)N3)[C@@H](C)O",
    "BGC0000543":r"CC[C@H](C)[C@@H](C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](C)C(=O)N[C@H](CS)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CS)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCCCN)C(=O)N/C(=C\C)/C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](C)C(=O)N/C(=C\C)/C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CC1=CC=CC=C1)C(=O)N/C(=C(\C)/S)/C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CS)C(=O)N[C@@H](CS)C(=O)N[C@@H](CCCCN)C(=O)NCC(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(=O)N)C(=O)NCC(=O)N[C@@H](CS)C(=O)N[C@@H](CCCCN)C(=O)O)NC(=O)[C@H](C)NC(=O)[C@@H]2CCCN2C(=O)CNC(=O)[C@H](C)NC(=O)C(=O)CC",
    "BGC0000529":r"O=C(N[C@H](C(NC(C(N[C@H](C(N[C@H](C(N[C@H]([C@H](C)SC[C@@H](C(N[C@H]([C@H](O)C)C(N[C@@H](C(N1[C@H]2[C@@H](O)CC1)=O)CSC[C@H](NC([C@H](CC(N)=O)NC([C@H](NC(CNC(CNC(CNC2=O)=O)=O)=O)CSC[C@@H](C(N/C=C\SC3)=O)NC4=O)=O)=O)C(N[C@@H]3C(N[C@@H]4CC5=CC=CC=C5)=O)=O)=O)=O)NC6=O)C(N7[C@H](C(NC6)=O)CCC7)=O)=O)CSC8)=O)CC(C)C)=O)=C)=O)CC9=CNC%10=C9C=C(Cl)C=C%10)[C@@H]8NC(/C(NC([C@H](C(C)C)N)=O)=C/C)=O",
    "BGC0000514":r"CC[C@H](C)[C@@H](C(=O)N[C@@H](C)C(=O)N[C@@H]1CSC[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC1=O)CCCCN)CC2=CC=CC=C2)CC(C)C)C(=O)N[C@@H]3[C@@H](SC[C@@H](NC(=O)CNC(=O)[C@@H]4CCCN4C3=O)C(=O)N[C@@H](C)C(=O)N[C@@H](CCCCN)C(=O)N/C(=C/C)/C(=O)NCC(=O)N[C@@H]5CSC[C@@H]6C(=O)N/C=C/SC[C@H](C(=O)N[C@H](C(=O)N6)CC7=CC=C(C=C7)O)NC(=O)[C@@H](NC(=O)[C@@H](NC5=O)CC8=CC=CC=C8)CC(=O)N)C)N",
    "BGC0000506":r"NC(CC1=CNC2=C1C=CC=C2)C(NC(CCCCN)C(NC(CSCC(C(NC(C5C)C(N4C(C(NC([H])C(NC(CS5)C(NC(C(C)C)C(NC(C6C)C(NC([H])C(NC(CC(C)C)C(NC(CC(C)C)C(NC(CCC(N)=O)C(NC(C(NC(CS6)C(NC(CC7=CC=CC=C7)C(NC(CC(C)C)C(NC(CCC(N)=O)C(NC(C8C)C(NC(C(CC)C)C(NC(C9C)C(NC(CS8)C(NC(CC(N)=O)C(NC(CS9)C(NC(CCCCN)C(NC(C(CC)C)C(NC(C(NC(CCCCN)C(O)=O)=O)=C)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=O)=CC)=O)=O)=O)=O)=O)=O)=O)=O)=O)CCC4)=O)=O)NC3=O)C(NC(CCC(O)=O)C(NC(C(NC3C(C)C)=O)=C)=O)=O)=O)=O",
    "BGC0002346":r"CC1CC(=O)C(CC2=CC=CC=C2)NC(=O)C2NC(=O)C(CC3=CC=CC=C3)NC(=O)C3CCCN3C(=O)CNC(=O)C(C)NC(=O)C3CNCCCCC(C(=O)O)NC(=O)C4NC(=O)C(CO)NC(=O)CNC(=O)C(C(O)C(=O)O)NC(=O)C(CSC(C)C(NC(=O)C(CO)NC(=O)C(C)NC(=O)C(N)CSC4C)C(=O)NC(CSC2C)C(=O)N3)NC1=O",
    "BGC0000539":r"CCC(C)C(C(=O)N1CCCC1C(=NC2C(SSCC3C(=NC(C(=NC(C(=NC(C(=NC4CSSCC(N=C(C(CSSCC(C(=NC(C(=NC(C(=N3)O)CC(=O)O)O)CC5=CN=CN5)O)N=C(C(N=C2O)C(C)C)O)N=C(C(=CC)N=C(C(N=C(C(N=C(C(N=C(C(N=C(C(N=C4O)CC6=CC=CC=C6)O)CCC(=N)O)O)CC7=CC=CC=C7)O)C(C)C)O)CC8=CC=CC=C8)O)O)O)C(=NC(CO)C(=O)O)O)O)CC(=N)O)O)CCSC)O)CC9=CN=CN9)O)C)O)N=C(C(C(C)C)N=C(CN=C(C(CO)N=C(C(CCCCN)N=C(C(CCCCN)N=C(C(CCCCN)N)O)O)O)O)O)O",
    "BGC0000527":r"CCC(C)C1C(=O)NC=CSC(C2C(=O)NC(=C)C(=O)NC(C(=O)NC(CSC(C(C(=O)NC(C(=O)N2)CC(C)C)NC(=O)C3CSC(C(C(=O)NC(C(=O)N4CCCC4C(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)NC(C(=O)N3)C(C)C)CC(C)C)NC(=O)C(CC5=CC=CC=C5)NC(=O)C6C(SCC(C(=O)N6)N)C)C)C)C(=O)N1)CCC(=O)O)C",
    "BGC0000495":r"CCC(C)C1C(=O)NC(C(=O)NC2CSC(C(C(=O)NC(C(=O)NC(C(SCC3C(=O)NC(C(=O)NC(CS(=O)C(C(C(=O)NC(C(=O)NC(C(=O)N3)C(C)CC)C(C)C)NC(=O)CNC2=O)C)C(=O)O)C)C)C(=O)N1)CC(C)C)NC(=O)C4CSCC(C(=O)NC(C(=O)NCC(=O)NC(C(=O)NC(C(=O)N4)C(C)C)CC5=CNC6=CC=CC=C65)CO)N)C)CCC(=O)O",
    "BGC0001860":r"O=C1C2=CSC(C3=C(C4=CSC(C5=CSC([C@@H](NC([C@H](CC6=CC=C(O)C=C6)NC(C7=CSC([C@H]([C@H](O)C8=CC=CC=C8)NC(C9=C(C)SC([C@H](CC(N)=O)N1)=N9)=O)=N7)=O)=O)[C@H](C)C%10CO%10)=N5)=N4)N=C(C%11=NC(C(NC(C(NC(C(N[C@@H](C)C(O)=O)=O)=C)=O)=C)=O)=CS%11)C=C3)=N2",
    "BGC0001474":r"C[C@H]([C@H](NC1=S)C2=NC([C@@H](N3)[C@](C4=NC(C(N[C@@](C(N/C(C5=N[C@@H](C(N[C@](C(C)(O)[C@@H](C)O)([H])C6=NC1=CS6)=O)CS5)=C\C)=O)([H])[C@@H](C)O)=O)=CS4)(NC([C@H](C)NC(C(NC([C@H](C)NC([C@H](C(C)C)N[C@@H]7C=CC8=C(N=C9C=C8[C@@H](O)C)[C@H]7O)=O)=O)=C)=O)=O)CC[C@@H]3C%10=NC(C(NC(C(NC(C(O)=O)=C)=O)=C)=O)=CS%10)=CS2)OC9=O",
    "BGC0001473":r"C[C@H]([C@H](NC1=S)C2=NC([C@@H](N3)[C@](C4=NC(C(N[C@@](C(N/C(C5=N[C@@H](C(N[C@](C(C)(O)[C@@H](C)O)([H])C6=NC1=CS6)=O)CS5)=C\C)=O)([H])[C@@H](C)O)=O)=CS4)(NC([C@H](C)NC(C(NC([C@H](C)NC([C@H]([C@@H](CC)C)N[C@@H]7C=CC8=C(N=C9C=C8[C@@H](O)C)[C@H]7O)=O)=O)=C)=O)=O)CC[C@@H]3C%10=NC(C(NC(C(NC(C(N)=O)=C)=O)=C)=O)=CS%10)=CS2)OC9=O",
    "BGC0001471":r"C[C@H]([C@H](NC1=S)C2=NC([C@](N3)([H])[C@](C4=NC(C(N[C@@](C(N/C(C5=N[C@](C(N[C@](C(C)(O)[C@@H](C)O)([H])C6=NC1=CS6)=O)([H])CS5)=C\C)=O)([H])[C@@H](C)O)=O)=CS4)(NC(C(NC(C(NC(C(NC([C@@](C)([H])N[C@@H]7C=CC8=C(N=C9C=C8[C@@H](O)C)[C@H]7O)=O)=C)=O)=C)=O)=C)=O)CCC3C%10=NC(C(NC(C(NC(C(O)=O)=C)=O)=C)=O)=CS%10)=CS2)OC9=O",
    "BGC0001155":r"CC1=C2C(=O)NCC(=O)NC(C3=NC(=CS3)C4=NC(=CS4)C5=C(C=CC(=N5)C6=NC(=CS6)C7=NC(CO7)C(=O)N8CCCC8C(=O)N)C9=NC(=CS9)C(=O)NC(C3=NC(=C(S3)COC)C(=O)NC(C(=N2)S1)CC(=O)NC)C(C)C)C(C1=CC=CC=C1)O",
    "BGC0001146":r"O=C([C@H]1CSC(C2=CSC(C3=CSC(C(NC([C@@H]4CSC(C5=CSC(C6=CC=C(C(N[C@@H](CC(N)=O)C7=N[C@H](C(N/C8=C\C)=O)CS7)=O)C([C@H](NC9=O)C)=N6)=N5)=N4)=O)=C)=N3)=N2)=N1)N[C@H](C(O)=O)CS[C@@](C)(C(N[C@@H](CCCCN)C(NCC(NC(C(N%10[C@H]9CCC%10)=O)=C)=O)=O)=O)NC8=O",
    "BGC0001145":r"[H][C@]12CSC(=N1)C1=CSC(=N1)C1=CSC(=N1)C(=C)NC(=O)[C@@]1([H])CSC(=N1)C1=CSC(=N1)C1=CC=C3C(=N1)[C@H](C)NC(=O)[C@]1([H])CCCN1C(=O)\C(NC(=O)CNC(=O)[C@H](CCCNC(N)=N)NC(=O)C(C)(NC(=O)\C(NC(=O)C1CSC(=N1)[C@H](CC(N)=O)NC3=O)=C\C)SC[C@H](NC2=O)C(O)=O)=C\C",
    "BGC0000628":r"O=C(NC(C(O)=O)C)C1=CSC(C2=CSC(C3=NC(C4=CSC(C(CSC(C5=CNC6=C5C=CC=C6)=O)NC7=O)=N4)=C(C8=NC(C(NC(C(C)O)C(N/C(C9=NC(C(NC(C%10=NC7CS%10)C(O)C)=O)=CS9)=C\C)=O)=O)=CS8)C=C3)=N2)=N1",
    "BGC0000615":r"C\C=C1\NC(=O)[C@@H](NC(=O)C2=CSC(=N2)C2=C(N=C(C=C2)C(=O)NC(=C)C(=O)NC(=C)C(=O)NCC(C)=O)C2=COC(=N2)C(=C)NC(=O)C(=C)NC(=O)C2=CSC(=N2)[C@H](C)NC(=O)C2=CSC(CNC(=O)C3=C(C)OC1=N3)=N2)[C@@H](C)O",
    "BGC0000614":r"CCC(C)C1C(=O)NC(C(=O)NC(=C)C(=O)NC(C(=O)NC23CCC(=NC2C4=CSC(=N4)C(C(OC(=O)C5=NC6=C(C=CC(C6O)N1)C(=C5)C(C)O)C)NC(=O)C7=CSC(=N7)C(NC(=O)C8CSC(=N8)/C(=C\C)/NC(=O)C(NC(=O)C9=CSC3=N9)C(C)O)C(C)(C(C)O)O)C1=NC(=CS1)C(=O)NC(=C)C(=O)NC(=C)C(=O)N)C)C",
    "BGC0000613":r"C[C@@H](C1CO1)[C@@H]1NC(=O)[C@H](CC2=CC=C(O)C=C2)NC(=O)C2=CSC(=N2)[C@@H](NC(=O)C2=C(C)SC(=N2)[C@H](CC(N)=O)NC(=O)C2=CSC(=N2)C2=C(N=C(C=C2)C2=NC(=CS2)C(=O)NC(=C)C(=O)NC(=C)C(O)=O)C2=CSC(=N2)C2=CSC1=N2)[C@H](O)C1=CC=CC=C1",
    "BGC0000612":r"C/C=C/1\C2=NC(=CS2)C(=O)NC(C3=NC(=CS3)C(=O)NC(C4=NC(=CS4)C5=C(C=CC(=N5)C6=NC(=CS6)C7=NC(=CS7)C(=O)N/C(=C/C)/C(=O)NCC(C)O)C8=NC(=CS8)C(=O)NC(C(=O)N1)C(C)O)C(C)O)C(C)(C)O",
    "BGC0000611":r"CC=C1C2=NC(CS2)C(=O)NC(C3=NC(=CS3)C(=O)NC4C(OC(=O)C5=NC6=C(C=CC(C6O)NC(C(=O)NC(=C)C(=O)NC(=C)C(=O)NC(C(=O)NC7(CCC(=NC7C8=CSC4=N8)C9=NC(=CS9)C(=O)NC(=C)C(=O)NC(=C)C(=O)N)C2=NC(=CS2)C(=O)NC(C(=O)N1)C(C)O)C)C(C)C)C(=C5)C(C)O)C)C(C)(C(C)O)O",
    "BGC0000610":r"C/C=C\1/C2=NC(=CS2)C(=O)NC3CC(C(=O)OCC4=C5C(=C(C(=O)SCC(C6=NC(=CS6)C7=N/C(=C/8\NC(=CS8)C(=O)NC(=C)C(=O)N)/C(=O)C=C7C9=NC(=CS9)C(=O)NC(C(=O)N1)C(C)O)NC(=O)C1=CSC3=N1)NC5=CC=C4)C)O",
    "BGC0000609":r"C[C@H]1[C@H]([C@@](C[C@@H](O1)O[C@H]2[C@@H]3[C@H]4C5=NC(=CS5)C(=O)N[C@@H](COC(=O)C6=C(CO3)C7=C(COC2=O)C=CC=C7N6O)C8=NC(=CS8)C9=N/C(=C/1\NC(=CS1)C(=O)NC(=C)C(=O)N)/C(=O)C=C9C1=NC(=CS1)C(=O)N[C@H](C(=O)N/C(=C(\C)/OC)/C1=NC(=CS1)C(=O)N4)[C@@H](C)O)(C)O)N(C)C.C(=O)(C(F)(F)F)O",
    "BGC0000608":r"C[C@H]1[C@H]([C@@](C[C@@H](O1)O[C@H]2[C@@H]3[C@H]4C5=NC(=CS5)C(=O)N[C@@H](COC(=O)C6=C(CO3)C7=C(COC2=O)C=CC=C7N6O)C8=NC(=CS8)C9=N/C(=C/1\NC(=CS1)C(=O)NC(=C)C(=O)N)/C(=O)C=C9C1=NC(=CS1)C(=O)N[C@H](C(=O)N/C(=C(\C)/OC)/C1=NC(=CS1)C(=O)N4)[C@@H](C)O)(C)O)N(C)C.C(=O)(C(F)(F)F)O",
    "BGC0000607":r"C/C=C/1\C2=NC(=CS2)C(=O)NC(C3=NC(=CS3)C(=O)NC(C4=NC(=CS4)C5=C(C=CC(=N5)C6=NC(=CS6)C7=NC(=CS7)C(=O)N/C(=C/C)/C(=O)NCC(C)O)C8=NC(=CS8)C(=O)NC(C(=O)N1)C(C)O)C(C)O)C(C)C",
    "BGC0000606":r"CC1C(=NC(=C)C2NC(=CO2)C3=C(C=CC(=N3)C4=NC(=CS4)C(=NC(C)C(=NC(CCC(=N)O)C(=O)N5CCCC5C(=O)O)O)O)C(=NC(C(=NCC(=NC(C6=NC(=CS6)C(=NC(=C)C7=NC(=CS7)C(=NC(C(=N1)O)CCC(=N)O)O)O)CO)O)O)CC8=CNC9=CC=CC=C98)O)O",
    "BGC0000605":r"CC1CC(O)N2C1C1=NC(=CS1)C1=NC(=CS1)C1=NC(=CC=C1C1=NC(=C(C)O1)C(=O)NC(CC(N)=O)C1=NC(=CS1)C(=O)NC(CC1=CC=CC=C1)C1=NC(CS1)C(=O)NC(CC1=CC=C(O)C=C1)C2=O)C1=NC(=CS1)C(=O)NC(=C)C(=O)NC(=C)C(O)=O",
    "BGC0000604":r"CC1=C2C(=O)NCC(=O)NC(C3=NC(=CS3)C4=NC(=CS4)C5=C(C=CC(=N5)C6=NC(=CS6)C7=NC(CO7)C(=O)N8CCCC8C(=O)N)C9=NC(=CS9)C(=O)NC(C3=NC(=C(S3)COC)C(=O)NC(C(=N2)S1)CC(=O)NC)C(C)C)C(C1=CC=CC=C1)O",
    "BGC0000603":r"[H][C@](C)(O)[C@]1([H])NC(=O)[C@@]2(C)NC(=O)\C(NC(=O)[C@]3([H])CSC(=N3)[C@]([H])(CC(N)=O)NC(=O)C3=CC=C(N=C3[C@@]([H])(C)NC(=O)[C@]3([H])CCCN3C(=O)\C(NC(=O)CNC1=O)=C\C)C1=NC(=CS1)C1=N[C@@]([H])(CS1)C(=O)NC(=C)C1=NC(=CS1)C1=NC(=CS1)C1=N[C@@]([H])(CS1)C(=O)N[C@@]([H])(CS2)C(O)=O)=C\C"
}

#Parameter for å avgjøre ringstørrelse
pRINGSIZE = 9
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

#Sjekker om det går an å gjøre C-term tailoring per LanD
def c_term_oxidative_decarboxylationable(node):
    temp_mol = node.molecule
    cys_c_term_query = Chem.MolFromSmarts("[NX3H1][C]([C][S])[C](=[O])[OX2H1]")
    Chem.SanitizeMol(temp_mol)
    return temp_mol.HasSubstructMatch(cys_c_term_query)

#Sjekker om det går an å gjøre heterosyklisering av Cys (thiopeptid)
def heterocyclizationable(node):
    temp_mol = node.molecule
    cys_query = Chem.MolFromSmarts("[C:1](=[O:2])[N:3][C:4]([C:5][SX2H:6])[C:7](=[O:8])")
    Chem.SanitizeMol(temp_mol)
    return temp_mol.HasSubstructMatch(cys_query)

def macrocyclizationable(node):
    temp_mol = node.molecule
    two_pi_query = Chem.MolFromSmarts("[N:1][C:2](=[CX3H2:3])")
    four_pi_query = Chem.MolFromSmarts("[C:9](=[CX3H2:10])([c:11]2[s:12][c:13][c:14][n:15]2)[N:16][C:17](=[O:18])[c:19]3[n:20][c:21][s:22][c:23]3")
    Chem.SanitizeMol(temp_mol)
    return temp_mol.HasSubstructMatch(two_pi_query) and temp_mol.HasSubstructMatch(four_pi_query)

def contains_thioenol_c_term(mol):
    temp_mol = mol
    thioenol_c_term_query = Chem.MolFromSmarts("[NX3H1:1][C:2]=[C:3][SX2H1:4]")
    Chem.SanitizeMol(temp_mol)
    return temp_mol.HasSubstructMatch(thioenol_c_term_query)

#NB bør testes i eget script. Funksjon som dehydrerer serine og threonine ved glutamylering
def glutamylation_elimination(node, model, strictness):
    if strictness == 0:
        serine_query = Chem.MolFromSmarts("[C:1]([CX4H2:2][OX2H:3])[C:4](=[O:5])")
        threonine_query = Chem.MolFromSmarts("[C:1]([C:2]([C:3])[O:4])[C:5](=[O:6])")
        serine_smarts = "[C:1]([CX4H2:2][OX2H:3])[C:4](=[O:5]) >> [CX4H0:1](=[CX3H2:2])[C:4](=[O:5]).[O:3]"
        threonine_smarts = "[C:1]([C:2]([C:3])[O:4])[C:5](=[O:6]) >> [C:1](=[C:2][C:3])[C:5](=[O:6]).[O:4]"
    elif strictness == 1:
        serine_query = Chem.MolFromSmarts("[C:1]([CX4H2:2][OX2H:3])[C,c:4]")
        threonine_query = Chem.MolFromSmarts("[C:1]([C:2]([C:3])[O:4])[C,c:5]")
        serine_smarts = "[C:1]([CX4H2:2][OX2H:3])[C,c:4] >> [CX4H0:1](=[CX3H2:2])[C,c:4].[O:3]"
        threonine_smarts = "[C:1]([C:2]([C:3])[O:4])[C,c:5] >> [C:1](=[C:2][C:3])[C,c:5].[O:4]"
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
        while len(products) > 1:
            Chem.SanitizeMol(products[0][0])
            products = ser_rxn.RunReactants([products[0][0], ])
            Chem.SanitizeMol(products[0][0])
            n += 1

        substrate = products[0][0]

    if substrate.HasSubstructMatch(threonine_query):
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

    m_node = Metabolite_node(dehydrated_peptide, node.name, node, dehydration_reaction, substrate)
    return m_node

#NB bør testes i eget script. Funksjon som dehydrerer serine og threonine ved fosforylering
def phosphorylation_dehydration(node, model, strictness):
    if strictness == 0:
        serine_query = Chem.MolFromSmarts("[C:1]([CX4H2:2][OX2H:3])[C:4](=[O:5])")
        threonine_query = Chem.MolFromSmarts("[C:1]([C:2]([C:3])[O:4])[C:5](=[O:6])")
        serine_smarts = "[C:1]([CX4H2:2][OX2H:3])[C:4](=[O:5]) >> [CX4H0:1](=[CX3H2:2])[C:4](=[O:5]).[O:3]"
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

        while len(products) > 1:
            Chem.SanitizeMol(products[0][0])
            products = ser_rxn.RunReactants([products[0][0], ])
            Chem.SanitizeMol(products[0][0])
            n += 1

        substrate = products[0][0]

    if substrate.HasSubstructMatch(threonine_query):
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
         model.metabolites.get_by_id("pi_c"): n, node.metabolite: -1.0, dehydrated_peptide: 1.0, model.metabolites.get_by_id("mg2_c"):-n,
         model.metabolites.get_by_id("mg2_c"):n})

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

#Only supports lanthionine/methyllanthionine formation, labionin/avionin not implemented (completely)
def cyclization(node, model, bgc_dict):
    substrate = node.molecule
    #Counts cycles (pun intended), each corresponds to creation of a water molecule
    n = 0
    product, c = recursive_cyclization(substrate, n, bgc_dict, 0)
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
        formula=lan_d_smiles,
        name=node.name + " LanD tailored peptide",
        compartment="c"
    )

    lan_d_tailoring_reaction = Reaction("lan_d_tailoring_reaction of " + node.name)
    lan_d_tailoring_reaction.name = "LanD tailoring reaction of " + node.name
    lan_d_tailoring_reaction.lower_bound = 0.
    lan_d_tailoring_reaction.upper_bound = 1000.
    lan_d_tailoring_reaction.add_metabolites(
        {node.metabolite: -1.0, lan_d_tailored_peptide: 1.0, model.metabolites.get_by_id("fmn_c"): -1.0, model.metabolites.get_by_id("fmnh2_c"): 1.0,
         model.metabolites.get_by_id("co2_c"): 1.0}
    )

    m_node = Metabolite_node(lan_d_tailored_peptide, node.name, node, lan_d_tailoring_reaction, substrate)
    return m_node

#Splitt opp i to reaksjoner for hver ko-faktor
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

    while len(products) > 1:
        Chem.SanitizeMol(products[0][0])
        products = cys_hc_rxn.RunReactants([products[0][0], ])
        Chem.SanitizeMol(products[0][0])
        n += 1

    substrate = products[0][0]

    hc_smiles = Chem.MolToSmiles(substrate)

    thioazole_peptide = Metabolite(
        id=node.name + "_thioazole_containing_peptide_c",
        formula=hc_smiles,
        name=node.name + " thioazole containing peptide intermediate",
        compartment="c"
    )
    thioazole_reaction = Reaction("thioazole heterocyclization of " + node.name)
    thioazole_reaction.name = "Thioazole heterocyclization reaction of " + node.name
    thioazole_reaction.lower_bound = 0.
    thioazole_reaction.upper_bound = 1000.
    thioazole_reaction.add_metabolites(
        {model.metabolites.get_by_id("atp_c"): -n, model.metabolites.get_by_id("h2o_c"):n, model.metabolites.get_by_id("adp_c"): n,
         model.metabolites.get_by_id("pi_c"): n, node.metabolite: -1.0, thioazole_peptide: 1.0, model.metabolites.get_by_id("fmn_c"): -n,
         model.metabolites.get_by_id("fmnh2_c"): n})

    m_node = Metabolite_node(thioazole_peptide, node.name, node, thioazole_reaction, substrate)
    return m_node

#Reactions that modify the thiopeptide central macrocycle by the type
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
            formula=smiles,
            name=node.name + " series e thiopeptide",
            compartment="c"
        )

        modification_reaction = Reaction("Series E modification of " + node.name)
        modification_reaction.name = "Series E modification reaction of " + node.name
        modification_reaction.lower_bound = 0.
        modification_reaction.upper_bound = 1000.
        modification_reaction.add_metabolites(
            {node.metabolite: -1.0, series_e_peptide: 1.0, model.metabolites.get_by_id("h2o_c"): -1.0}
        )
        m_node = Metabolite_node(series_e_peptide, node.name, node, modification_reaction, product)
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
            formula=smiles,
            name=node.name + " series b thiopeptide",
            compartment="c"
        )

        modification_reaction = Reaction("Series B modification of " + node.name)
        modification_reaction.name = "Series B modification reaction of " + node.name
        modification_reaction.lower_bound = 0.
        modification_reaction.upper_bound = 1000.
        modification_reaction.add_metabolites(
            {node.metabolite: -1.0, series_b_peptide: 1.0, model.metabolites.get_by_id("h_c"): -1.0}
        )
        m_node = Metabolite_node(series_b_peptide, node.name, node, modification_reaction, product)
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
            formula=smiles,
            name=node.name + " series d thiopeptide",
            compartment="c"
        )

        modification_reaction = Reaction("Series D modification of " + node.name)
        modification_reaction.name = "Series D modification reaction of " + node.name
        modification_reaction.lower_bound = 0.
        modification_reaction.upper_bound = 1000.
        modification_reaction.add_metabolites(
            {node.metabolite: -1.0, series_d_peptide: 1.0, model.metabolites.get_by_id("h_c"): 1.0}
        )
        m_node = Metabolite_node(series_d_peptide, node.name, node, modification_reaction, product)
        return m_node
    else:
        return node

#Bruk ca. samme metodikk som i lanthipeptid sykliseringen, men trenger ikke være rekursivt (bare én ring)
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
        formula=cyc_smiles,
        name=node.name + " macrocyclic thiopeptide",
        compartment="c"
    )

    cyclization_reaction = Reaction("macrocyclization of " + node.name)
    cyclization_reaction.name = "Macrocyclization reaction of " + node.name
    cyclization_reaction.lower_bound = 0.
    cyclization_reaction.upper_bound = 1000.
    cyclization_reaction.add_metabolites(
        {node.metabolite: -1.0, cyclic_peptide: 1.0, model.metabolites.get_by_id("h2o_c"): 1}
    )

    m_node1 = Metabolite_node(cyclic_peptide, node.name, node, cyclization_reaction, cyc_mol)

    m_node2 = macrocycle_modification(m_node1, model, bgc_dict)

    return m_node2

def match(predicted_smiles, BGC_ID):
    mols = [Chem.MolFromSmiles(predicted_smiles), Chem.MolFromSmiles(BGC_to_smiles[BGC_ID])]
    fps = [Chem.RDKFingerprint(x) for x in mols]
    return DataStructs.FingerprintSimilarity(fps[0], fps[1])

def get_gb_list_from_antismash_output(cluster_path):  # yes
    # domains, pfam_entries, smCOG, EC_number, rxn_in, rxn_out, core_gene, amino_acid_sequence, predicted_EC, length
    # triple pound signs means that information is stored within a custom class object.
    gb_list = []
    for gb_record in SeqIO.parse(open(cluster_path, "r"), "genbank"):
        gb_list.append(gb_record)
    return gb_list

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

def run(bgc_path, output, mode):
    print("Running...")
    summary = pd.DataFrame({"Name":[], "Core Name":[], "Core Class":[], "Core Subclass":[], "Core Number":[], "BGC":[], "Organism":[], "End product SMILES":[], "Metabolite list":[], "Prediction Score":[]})
    column_list = ["Name", "Core Name", "Core Class", "Core Subclass", "Core Number", "BGC", "Organism", "End product SMILES", "Metabolite list", "Prediction Score"]
    bgc_path = Path(bgc_path)
    result_list = []

    if bgc_path.is_dir():
        for filename in os.listdir(bgc_path):
            p = os.path.join(bgc_path, filename)
            if os.path.isfile(p):
                print(filename)
                bgc_dict = {
                "Name":os.path.splitext(filename)[0], "Core Name":None, "Core Number":None, "BGC":"Unknown", "Organism":"Unknown", "End product SMILES":"Unknown", "Metabolite list":[],
                "SMCOG1070":False, "SMCOG1140":False, "SMCOG1155":False, "SMCOG1075":False, "SMCOG1030":False, "Core":None, "Leader":None, "LanO":False, "Prediction Score":None,
                "LanD":False, "Core Class":"Unknown", "Core Subclass":"Unknown", "YcaO":False, "Oxazole":False
                }
                result_list.extend(analyse_bgc(p, bgc_dict, output, mode))
    else:
        name = os.path.basename(bgc_path)
        bgc_dict = {
            "Name": os.path.splitext(name)[0], "Core Number":None, "Core Name":None, "BGC": "Unknown", "Organism": "Unknown", "End product SMILES": "Unknown", "Metabolite list": [],
            "SMCOG1070": False, "SMCOG1140": False, "SMCOG1155": False, "SMCOG1075": False, "SMCOG1030":False, "Core": None, "Leader": None, "LanO":False, "Prediction Score":None,
            "LanD":False, "Core Class":"Unknown", "Core Subclass":"Unknown", "YcaO":False, "Oxazole":False
        }
        result_list.extend(analyse_bgc(bgc_path, bgc_dict, output, mode))

    for result in result_list:
        new_row = {}
        for column in column_list:
            new_row[column] = [result[column]]
        temp_df = pd.DataFrame(new_row)
        summary = pd.concat([temp_df, summary.loc[:]]).reset_index(drop=True)

    summary.to_csv(os.path.join(output, "thiopeptide performance.csv"))

def parse_gbk(bgc_path, bgc_dict):
    count = 0
    gb_list = get_gb_list_from_antismash_output(bgc_path)
    for gb_record in gb_list:
        bgc_dict["Organism"] = gb_record.annotations["organism"]
        for feat in gb_record.features:
            if feat.type == "CDS":
                if feat.qualifiers.get("gene_functions") != None:
                    for i in feat.qualifiers.get("gene_functions"):
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
    bgc_dict = add_to_model(bgc_dict, output, mode)
    return bgc_dict

def add_to_model(bgc_dict, output, mode):
    model = cobra.io.read_sbml_model(
        join("C:\\Users\\adtho\\Documents\\Python Scripts\\pythonProject-20221219T143334Z-001\\pythonProject\\models",
             "SCO-GEM.xml"))

    leaderseq = bgc_dict["Leader"][0]
    coreseq = bgc_dict["Core"][0]
    pepseq = leaderseq+coreseq

    origin_node = translation(pepseq, bgc_dict["Name"], model)
    current_node = origin_node
    current_node = core_cleave(current_node, coreseq, leaderseq, model)
    if heterocyclizationable(current_node) and bgc_dict["YcaO"] or heterocyclizationable(current_node) and bgc_dict["Core Class"] == "thiopeptide":
        current_node = heterocyclization(current_node, model, bgc_dict["Core Subclass"])
    #Glutamylation dehydration of either class I lanthipeptides or thiopeptides
    if dehydratable(current_node) and bgc_dict["SMCOG1155"] or dehydratable(current_node) and bgc_dict["Core Class"] == "thiopeptide":
        current_node = glutamylation_elimination(current_node, model, 1)
    #Phosphorylation dehdyration of class II, III or IV lanthipeptides
    if dehydratable(current_node) and bgc_dict["SMCOG1030"] or dehydratable(current_node) and bgc_dict["SMCOG1070"]:
        current_node = phosphorylation_dehydration(current_node, model, 0)
    #LanD mediated C-terminus tailoring of class I lanthipeptides
    if c_term_oxidative_decarboxylationable(current_node) and bgc_dict["LanD"] and bgc_dict["Core Class"] == "lanthipeptide":
        current_node = lan_d_tailoring(current_node, model)
    #Cyclization of lanthipeptides
    #if cyclizationable(current_node) and bgc_dict["SMCOG1140"] or cyclizationable(current_node) and bgc_dict["SMCOG1070"] or cyclizationable(current_node) and bgc_dict["SMCOG1030"]:
    if cyclizationable(current_node):
        current_node = cyclization(current_node, model, bgc_dict)
    #LanO mediated N-terminus tailoring of class I lanthipeptides
    if n_term_lac_able(current_node) and bgc_dict["LanO"] and bgc_dict["Core Class"] == "lanthipeptide":
        current_node = lan_o_tailoring(current_node, model)
    # Primary macrocyclization of thiopeptide (Need to implement secondary macrocycle or siderings)
    if macrocyclizationable(current_node) and bgc_dict["Core Class"] == "thiopeptide":
        current_node = macrocyclization(current_node, model, bgc_dict)

    bgc_dict["End product SMILES"] = current_node.metabolite.formula

    if mode == 1 and bgc_dict["Name"] in BGC_to_smiles.keys():
        bgc_dict["Prediction Score"] = match(bgc_dict["End product SMILES"], bgc_dict["Name"])
        print(bgc_dict["Prediction Score"])

    metabolite_list = []

    while current_node != None:
        print(current_node.reaction)
        metabolite_list.append(current_node.metabolite.id)
        current_node = current_node.parent
    metabolite_list.reverse()

    bgc_dict["Metabolite list"] = metabolite_list
    bgc_dict["Success"] = "Yes"
    #Missing functionality: Output sbml file for model, some sort of image handling with appropriate naming, fingerprint molecular matching for characterized products
    return bgc_dict

def main():
    run("C:\\Users\\adtho\\Documents\\Python Scripts\\pythonProject-20221219T143334Z-001\\pythonProject\\antiSMASH 6 annotated RiPP clusters (lanthipeptides+thiopeptides)",
        "C:\\Users\\adtho\\Documents\\Python Scripts\\pythonProject-20221219T143334Z-001\\pythonProject\\Output", 1)

if __name__ == "__main__":
    main()