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

def match(generated_smiles, true_smiles, name):
    mols = [Chem.MolFromSmiles(generated_smiles), Chem.MolFromSmiles(true_smiles)]
    fps = [Chem.RDKFingerprint(x) for x in mols]
    print("Molecular similarity between generated "+name+" and characterized "+name+" is",
          DataStructs.FingerprintSimilarity(fps[0], fps[1]))

#Test for epilancin 15x, with and without tailoring
print("With tailoring")
match(r"[CH2][C](O)C(=O)[NH:6][C@H:7]([C:8](=[O:9])[NH:11]C(=C)C(=O)[NH:17][C@H:18]([C:19](=[O:20])[NH:25][C@H:26]([C:27](=[O:28])[NH:32][C@H:33]([C:34](=[O:35])[NH:41]C(=CC)C(=O)[NH:48]C1C(=O)[NH:55][C@@H:56]([C@@H:59]([CH3:60])[CH2:61][CH3:62])[C:57](=[O:58])[NH:63][C@@H:64]([CH2:67][CH2:68][CH2:69][CH2:70][NH2:71])[C:65](=[O:66])[NH:72][C@@H:73]([CH3:76])[C:74](=[O:75])[NH:77][C]2CSC[C@@H](C(=O)[NH:115][C@@H:116]([CH2:119][CH2:120][CH2:121][NH:122][C:123](=[NH:124])[NH2:125])[C:117](=[O:118])[NH:126][CH2:127][C:128](=[O:129])[NH:130][C@@H:131]([CH2:134][c:135]3[cH:136][cH:138][cH:140][cH:139][cH:137]3)[C:132](=[O:133])[NH:141]C3C(=O)[NH:148][C@@H:149]([CH2:152][CH:153]([CH3:154])[CH3:155])[C:150](=[O:151])[NH:156]C(=CC)C(=O)[NH:163][C@@H](CSC1C)C(=O)[NH:169][CH2:170][C:171](=[O:172])[NH:173][C@H](C(=O)[NH:179][C@H:180]([C:181](=[O:182])[NH:189][C@H:190]([C:191](=[O:192])[NH:200]C(=CC)C(=O)[NH:207][CH2:208][C:209](=[O:210])[NH:211][C@H:212]([C:213](=[O:214])[NH:220][C@H:221]([C:222](=[O:223])[OH:229])[CH2:224][CH2:225][CH2:226][CH2:227][NH2:228])[CH2:215][CH2:216][CH2:217][CH2:218][NH2:219])[CH2:193][c:194]1[cH:195][cH:197][cH:199][cH:198][cH:196]1)[CH2:183][c:184]1[n:185][cH:187][nH:188][cH:186]1)CSC3C)[NH:109][C:103](=[O:104])[C@H:102]([CH2:105][CH:106]([CH3:107])[CH3:108])[NH:101][C:94](=[O:95])[C@H:93]([CH2:96][CH2:97][CH2:98][CH2:99][NH2:100])[NH:92][C:85](=[O:86])[C@H:84]([CH2:87][CH2:88][CH2:89][CH2:90][NH2:91])[NH:83]C2=O)[CH2:36][CH2:37][CH2:38][CH2:39][NH2:40])[CH:29]([CH3:30])[CH3:31])[C@@H:21]([CH3:22])[CH2:23][CH3:24])[CH3:10]",
      r"CC[C@H](C)[C@H](NC(=O)C(\NC(=O)C(\NC(=O)[C@H](CCCCN)NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)C(=C)NC(=O)[C@H](C)NC(=O)C(C)O)[C@@H](C)CC)C(C)C)=C/C)=C/C)C(=O)NC(CCCCN)C(=O)N[C@@H](C)C(=O)N[C@H]1CSC[C@H](NC(=O)[C@H](CC(C)C)NC(=O)C(CCCCN)NC(=O)C(CCCCN)NC1=O)C(=O)N[C@@H](CCCNC(N)=N)C(=O)NCC(=O)N[C@@H](CC1=CC=CC=C1)C(=O)N[C@H]1C(C)SC[C@@H]2NC(=O)[C@@H](NC(=O)[C@H](CC(C)C)NC1=O)C(C)SC[C@H](NC(=O)CNC2=O)C(=O)N[C@@H](CC1=CN=CN1)C(=O)N[C@@H](CC1=CC=CC=C1)C(=O)N\C(=C\C)C(=O)NCC(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(O)=O",
      "Epilancin 15x")
print("Without tailoring")
match(r"C=C(N)C(=O)[NH:6][C@H:7]([C:8](=[O:9])[NH:11]C(=C)C(=O)[NH:17][C@H:18]([C:19](=[O:20])[NH:25][C@H:26]([C:27](=[O:28])[NH:32][C@H:33]([C:34](=[O:35])[NH:41]C(=CC)C(=O)[NH:48]C1C(=O)[NH:55][C@@H:56]([C@@H:59]([CH3:60])[CH2:61][CH3:62])[C:57](=[O:58])[NH:63][C@@H:64]([CH2:67][CH2:68][CH2:69][CH2:70][NH2:71])[C:65](=[O:66])[NH:72][C@@H:73]([CH3:76])[C:74](=[O:75])[NH:77][C]2CSC[C@@H](C(=O)[NH:115][C@@H:116]([CH2:119][CH2:120][CH2:121][NH:122][C:123](=[NH:124])[NH2:125])[C:117](=[O:118])[NH:126][CH2:127][C:128](=[O:129])[NH:130][C@@H:131]([CH2:134][c:135]3[cH:136][cH:138][cH:140][cH:139][cH:137]3)[C:132](=[O:133])[NH:141]C3C(=O)[NH:148][C@@H:149]([CH2:152][CH:153]([CH3:154])[CH3:155])[C:150](=[O:151])[NH:156]C(=CC)C(=O)[NH:163][C@@H](CSC1C)C(=O)[NH:169][CH2:170][C:171](=[O:172])[NH:173][C@H](C(=O)[NH:179][C@H:180]([C:181](=[O:182])[NH:189][C@H:190]([C:191](=[O:192])[NH:200]C(=CC)C(=O)[NH:207][CH2:208][C:209](=[O:210])[NH:211][C@H:212]([C:213](=[O:214])[NH:220][C@H:221]([C:222](=[O:223])[OH:229])[CH2:224][CH2:225][CH2:226][CH2:227][NH2:228])[CH2:215][CH2:216][CH2:217][CH2:218][NH2:219])[CH2:193][c:194]1[cH:195][cH:197][cH:199][cH:198][cH:196]1)[CH2:183][c:184]1[n:185][cH:187][nH:188][cH:186]1)CSC3C)[NH:109][C:103](=[O:104])[C@H:102]([CH2:105][CH:106]([CH3:107])[CH3:108])[NH:101][C:94](=[O:95])[C@H:93]([CH2:96][CH2:97][CH2:98][CH2:99][NH2:100])[NH:92][C:85](=[O:86])[C@H:84]([CH2:87][CH2:88][CH2:89][CH2:90][NH2:91])[NH:83]C2=O)[CH2:36][CH2:37][CH2:38][CH2:39][NH2:40])[CH:29]([CH3:30])[CH3:31])[C@@H:21]([CH3:22])[CH2:23][CH3:24])[CH3:10]",
      r"CC[C@H](C)[C@H](NC(=O)C(\NC(=O)C(\NC(=O)[C@H](CCCCN)NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)C(=C)NC(=O)[C@H](C)NC(=O)C(C)O)[C@@H](C)CC)C(C)C)=C/C)=C/C)C(=O)NC(CCCCN)C(=O)N[C@@H](C)C(=O)N[C@H]1CSC[C@H](NC(=O)[C@H](CC(C)C)NC(=O)C(CCCCN)NC(=O)C(CCCCN)NC1=O)C(=O)N[C@@H](CCCNC(N)=N)C(=O)NCC(=O)N[C@@H](CC1=CC=CC=C1)C(=O)N[C@H]1C(C)SC[C@@H]2NC(=O)[C@@H](NC(=O)[C@H](CC(C)C)NC1=O)C(C)SC[C@H](NC(=O)CNC2=O)C(=O)N[C@@H](CC1=CN=CN1)C(=O)N[C@@H](CC1=CC=CC=C1)C(=O)N\C(=C\C)C(=O)NCC(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(O)=O",
      "Epilancin 15x")
#print("RiPPminer Genome")
#Nettsiden deres er foreløpig nede