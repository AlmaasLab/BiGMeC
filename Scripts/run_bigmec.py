from Bio import SeqIO
import re
import os
import json
import test
import csv
import sys
import re
import pandas as pd
import cobra  
import plotly.graph_objects as go
import copy
from itertools import groupby
import warnings


'''
Change the 5 following paths:
'''
biggbk = "/Users/fredrikfossheim/Desktop/Master/master/antismash_output/webscraped_fixed/"
# 1) Folder containing all gbk files you want to translate into metabolic pathways
#    They are saved as models, and can be merged with the GEM later.
#    In this repository, the models that are found in "gbk_db_output_models.zip" are the pathways that 
#    Have been constructed. In that case, this folder contained all antiSMASH output files found in 
#    antismash_output_mibig_gbk_files1.zip
#    antismash_output_mibig_gbk_files2.zip
#    and
#    antismash_output_mibig_gbk_files3.zip
output_gbk = "/Users/fredrikfossheim/Desktop/Master/master/gbk_db_output_models/"
# 2) Folder that regular models are output (empty folder)
output_gbk_lump = "/Users/fredrikfossheim/Desktop/Master/master/gbk_db_output_models_lump/"
# 3) Folder that lump models are output (empty folder)
json_folder = '/Users/fredrikfossheim/Desktop/Master/master/json_files/difficidin/'
# 4) Folder that json files are output (empty folder)(an unnecessary step, but it may help look at the information that is saved)
sco_metabolic_model_path = '/Users/fredrikfossheim/Desktop/Master/master/streptomyces_coelicolor_metabolic_model.xml'
# 5) path of the genome scale metabolic model
snorre_model = cobra.io.read_sbml_model(sco_metabolic_model_path)

reducing_domains = [
    "PKS_DH",
    "PKS_DH2",
    "PKS_DHt",
    "PKS_KR",
    "PKS_ER",
    "MT",
    "Thioesterase"
]

dh_er_domains = [
    "PKS_DH",
    "PKS_DH2",
    "PKS_DHt",
    "PKS_ER"
]

mt_domains = [
    'cMT',
    'oMT',
    'nMT'
]

loader_domains = [
    'CAL_domain',
    'FkbH',
    'GNAT'
]

loader_at_domains = [
    'AMP-binding',
    'PKS_AT'
]

loader_acp_domains = [
    "ACP",
    "PP-binding",
    'ACP_beta',
    'PCP'
]

exclude_modules = [
    'DHD',
    'oMT'
]

general_domain_dict = {
    "PKS_DH": "PKS_DH",
    "PKS_DH2": "PKS_DH",
    "PKS_DHt": "PKS_DH",
    "PKS_KR": "PKS_KR",
    "PKS_ER": "PKS_ER",
    "MT": "PKS_cMT",
    "cMT": "cMT",
    "oMT": "oMT",
    "nMT": "nMT"
}

acid_to_bigg = {'gly': 'gly_c',
                'ala': 'ala__L_c',
                'val': 'val__L_c',
                'leu': 'leu__L_c',
                'ile': 'ile__L_c',
                'asp': 'asp__L_c',
                'asn': 'asn__L_c',
                'glu': 'glu__L_c',
                'gln': 'gln__L_c',
                'ser': 'ser__L_c',
                'thr': 'thr__L_c',
                'met': 'met__L_c',
                'cys': 'cys__L_c',
                'lys': 'lys__L_c',
                'arg': 'arg__L_c',
                'his': 'his__L_c',
                'pro': 'pro__L_c',
                'phe': 'phe__L_c',
                'tyr': 'tyr__L_c',
                'trp': 'trp__L_c'}

cofactor_metabolites_dict = {
    'nadph': snorre_model.metabolites.get_by_id('nadph_c'),
    'nadp': snorre_model.metabolites.get_by_id('nadp_c'),
    'h2o': snorre_model.metabolites.get_by_id('h2o_c'),
    'coa': snorre_model.metabolites.get_by_id('coa_c'),
    'co2': snorre_model.metabolites.get_by_id('co2_c'),
    'h+': snorre_model.metabolites.get_by_id('h_c'),
    'sam': snorre_model.metabolites.get_by_id('amet_c'),
    'sah': snorre_model.metabolites.get_by_id('ahcys_c'),
    'mxcoa': cobra.Metabolite('mxmal_c', formula='C14H20N6O5S', name='Methoxymalonyl-CoA', compartment='c'),
    '4hbf': cobra.Metabolite('4hbf_c', formula='na', name='4-hydroxy-benzoyl-formate', compartment='c'),
    'hpg': cobra.Metabolite('4hpg_c', formula='na', name='4-hydroxy-phenyl-glycine', compartment='c'),
    'dpg': cobra.Metabolite('dpg_c', formula='na', name='dihydroxy-phenyl-glycine', compartment='c'),
    'bht': cobra.Metabolite('bht_c', formula='na', name='beta-hydroxy-tyrosine', compartment='c'),
    'pip': cobra.Metabolite('Lpipecol_c', formula='na', name='pipecolic acid', compartment='c')
}

cofactor_reactions_dict = {  # reactions that are specific to certain domains
    'PKS_KR': {cofactor_metabolites_dict['nadph']: -1.0, cofactor_metabolites_dict['h+']: -1.0,
               cofactor_metabolites_dict['nadp']: 1.0},
    'cMT': {cofactor_metabolites_dict['sam']: -1.0, cofactor_metabolites_dict['sah']: 1.0},
    'oMT': {cofactor_metabolites_dict['sam']: -1.0, cofactor_metabolites_dict['sah']: 1.0},
    'PKS_DH': {cofactor_metabolites_dict['h2o']: 1.0},
    'PKS_ER': {cofactor_metabolites_dict['nadph']: -1.0, cofactor_metabolites_dict['h+']: -1.0,
               cofactor_metabolites_dict['nadp']: 1.0},
    'PKS_TE': {cofactor_metabolites_dict['h2o']: -1.0},
    'PKS_AT': {cofactor_metabolites_dict['coa']: 1.0, cofactor_metabolites_dict['co2']: 1.0},
    'Condensation': {cofactor_metabolites_dict['h2o']: -1.0},
    'nMT': {cofactor_metabolites_dict['sam']: -1.0, cofactor_metabolites_dict['sah']: 1.0},
    'mxmal': {snorre_model.metabolites.get_by_id('coa_c'): -1.0,
              snorre_model.metabolites.get_by_id('13dpg_c'): -1.0,
              snorre_model.metabolites.get_by_id('pi_c'): 2.0,
              snorre_model.metabolites.get_by_id('nadp_c'): -1.0,
              snorre_model.metabolites.get_by_id('nadph_c'): 1.0,
              snorre_model.metabolites.get_by_id('h_c'): 1.0,
              snorre_model.metabolites.get_by_id('amet_c'): -1.0,
              snorre_model.metabolites.get_by_id('ahcys_c'): 1.0,
              snorre_model.metabolites.get_by_id('fad_c'): -1.0,
              snorre_model.metabolites.get_by_id('fadh2_c'): 1.0,
              cofactor_metabolites_dict['mxcoa']: 1.0},

    'hpg_1': {snorre_model.metabolites.get_by_id('pphn_c'): -1,
              snorre_model.metabolites.get_by_id('34hpp_c'): 1,
              snorre_model.metabolites.get_by_id('co2_c'): 1,
              snorre_model.metabolites.get_by_id('h2o_c'): 1},

    'hpg_2': {snorre_model.metabolites.get_by_id('34hpp_c'): -1,
              snorre_model.metabolites.get_by_id('o2_c'): -1,
              snorre_model.metabolites.get_by_id('h2o_c'): 1,
              snorre_model.metabolites.get_by_id('4hmda_c'): 1},

    'hpg_3': {snorre_model.metabolites.get_by_id('4hmda_c'): -1,
              snorre_model.metabolites.get_by_id('fmn_c'): -1,
              snorre_model.metabolites.get_by_id('fmnh2_c'): 1,
              snorre_model.metabolites.get_by_id('nadh_c'): -1,
              snorre_model.metabolites.get_by_id('nad_c'): 1,
              cofactor_metabolites_dict['4hbf']: 1},

    'hpg_4': {cofactor_metabolites_dict['4hbf']: -1,
              snorre_model.metabolites.get_by_id('tyr__L_c'): -1,
              snorre_model.metabolites.get_by_id('34hpp_c'): 1,
              cofactor_metabolites_dict['hpg']: 1},

    'bht': {snorre_model.metabolites.get_by_id('tyr__L_c'): -1,
            snorre_model.metabolites.get_by_id('o2_c'): -1,
            snorre_model.metabolites.get_by_id('nadph_c'): -1,
            snorre_model.metabolites.get_by_id('h_c'): -1,
            snorre_model.metabolites.get_by_id('nadp_c'): 1,
            cofactor_metabolites_dict['hpg']: 1},

    'dpg': {snorre_model.metabolites.get_by_id('accoa_c'): -1,
            snorre_model.metabolites.get_by_id('malcoa_c'): -3,
            snorre_model.metabolites.get_by_id('coa_c'): 4,
            snorre_model.metabolites.get_by_id('co2_c'): 3,
            snorre_model.metabolites.get_by_id('h2o_c'): 1,
            snorre_model.metabolites.get_by_id('tyr__L_c'): -1,
            snorre_model.metabolites.get_by_id('34hpp_c'): 1,
            cofactor_metabolites_dict['dpg']: 1},

    'gnat': {snorre_model.metabolites.get_by_id('malcoa_c'): -1,
             snorre_model.metabolites.get_by_id('co2_c'): 1,
             snorre_model.metabolites.get_by_id('coa_c'): 1},

    'fkbh': {snorre_model.metabolites.get_by_id('13dpg_c'): -1.0,
             snorre_model.metabolites.get_by_id('pi_c'): 2.0},

    'ahba': {snorre_model.metabolites.get_by_id('e4p_c'): -1.0,
             snorre_model.metabolites.get_by_id('pep_c'): -1.0,
             snorre_model.metabolites.get_by_id('pi_c'): 1.0,
             snorre_model.metabolites.get_by_id('h2o_c'): 2.0},

    'acetyl': {snorre_model.metabolites.get_by_id('accoa_c'): -1,
               snorre_model.metabolites.get_by_id('coa_c'): 1,
               snorre_model.metabolites.get_by_id('co2_c'): 1},

    'shikimic_acid': {snorre_model.metabolites.get_by_id('skm_c'): -1,
                      snorre_model.metabolites.get_by_id('co2_c'): 1},

    'fatty_acid': {snorre_model.metabolites.get_by_id('accoa_c'): -1,
                   snorre_model.metabolites.get_by_id('malcoa_c'): -3,
                   snorre_model.metabolites.get_by_id('co2_c'): 4,
                   snorre_model.metabolites.get_by_id('coa_c'): 4},

    'NH2': {snorre_model.metabolites.get_by_id('malcoa_c'): -1,
            snorre_model.metabolites.get_by_id('co2_c'): 2,
            snorre_model.metabolites.get_by_id('gly_c'): -1,
            snorre_model.metabolites.get_by_id('coa_c'): 1
            },

    'pip': {snorre_model.metabolites.get_by_id('pyr_c'): -1,
            snorre_model.metabolites.get_by_id('lys__L_c'): -1,
            snorre_model.metabolites.get_by_id('h2o_c'): 1,
            snorre_model.metabolites.get_by_id('nadph_c'): -1,
            snorre_model.metabolites.get_by_id('nadp_c'): 1,
            snorre_model.metabolites.get_by_id('h_c'): -1,
            cofactor_metabolites_dict['pip']: 1}


}

tailoring_metabolites_dict = {
    'glucose_6_phosphate': snorre_model.metabolites.get_by_id('g6p_c'),
    '13biphosphoglycerate': snorre_model.metabolites.get_by_id('13dpg_c'),
    'succinyl_coa': snorre_model.metabolites.get_by_id('succoa_c'),
    'glycine': snorre_model.metabolites.get_by_id('gly_c'),
    'coa': snorre_model.metabolites.get_by_id('coa_c'),
    'co2': snorre_model.metabolites.get_by_id('co2_c'),
    'atp': snorre_model.metabolites.get_by_id('atp_c'),
    'amp': snorre_model.metabolites.get_by_id('amp_c'),
    'ppi': snorre_model.metabolites.get_by_id('ppi_c'),
    'pi': snorre_model.metabolites.get_by_id('pi_c'),
    'nadh': snorre_model.metabolites.get_by_id('nadph_c'),
    'nad': snorre_model.metabolites.get_by_id('nadp_c'),
    'h+': snorre_model.metabolites.get_by_id('h_c'),
    'h2o': snorre_model.metabolites.get_by_id('h2o_c')
}

tailoring_reactions_dict = {  # remember to add the secondary metabolite to these reactions
    'glycosyltransferase': {tailoring_metabolites_dict['glucose_6_phosphate']: -1,
                            tailoring_metabolites_dict['h+']: 1,
                            tailoring_metabolites_dict['pi']: 1
                            },
    'glycerol': {tailoring_metabolites_dict['13biphosphoglycerate']: -1.0,
                 tailoring_metabolites_dict['pi']: 2.0,
                 tailoring_metabolites_dict['h+']: -1.0,
                 tailoring_metabolites_dict['nad']: 1.0,
                 tailoring_metabolites_dict['nadh']: -1.0,
                 },
    'ALA': {tailoring_metabolites_dict['succinyl_coa']: -1.0,
            tailoring_metabolites_dict['glycine']: -1.0,
            tailoring_metabolites_dict['atp']: -1.0,
            tailoring_metabolites_dict['ppi']: 1.0,
            tailoring_metabolites_dict['co2']: 1.0,
            tailoring_metabolites_dict['coa']: 1.0,
            tailoring_metabolites_dict['amp']: 1.0,
            tailoring_metabolites_dict['h2o']: 1.0
            }

}

long_to_short = {'Malonyl-CoA': 'mal', 'Methylmalonyl-CoA': 'mmal', 'Methoxymalonyl-CoA': 'mxmal',
                 'Ethylmalonyl-CoA': 'emal', 'Isobutyryl-CoA': 'isobut', '2-Methylbutyryl-CoA': '2metbut',
                 'trans-1,2-CPDA': 'trans-1,2-CPDA', 'Acetyl-CoA': 'Acetyl-CoA', 'Benzoyl-CoA': 'benz',
                 'Propionyl-CoA': 'prop', '3-Methylbutyryl-CoA': '3metbut',
                 'CE-Malonyl-CoA': 'cemal', '2-Rhyd-Malonyl-CoA': '2Rhydmal', 'CHC-CoA': 'CHC-CoA',
                 'inactive': 'inactive'}

long_to_bigg = {'Malonyl-CoA': 'malcoa_c',
                'Methylmalonyl-CoA': 'mmcoa__R_c',
                'Methoxymalonyl-CoA': 'mxmal_c',
                'Ethylmalonyl-CoA': 'emcoa__S_c',
                'Isobutyryl-CoA': 'ibcoa_c',
                '2-Methylbutyryl-CoA': '2metbut',
                'trans-1,2-CPDA': 'trans-1,2-CPDA',
                'Acetyl-CoA': 'accoa',
                'Benzoyl-CoA': 'benzcoa',
                'Propionyl-CoA': 'ppcoa_c',
                '3-Methylbutyryl-CoA': '3metbut',
                'CE-Malonyl-CoA': 'cemal',
                '2-Rhyd-Malonyl-CoA': '2Rhydmal',
                'CHC-CoA': 'CHC-CoA',
                'inactive': 'inactive',
                'gly': 'gly_c',
                'ala': 'ala__L_c',
                'val': 'val__L_c',
                'leu': 'leu__L_c',
                'ile': 'ile__L_c',
                'asp': 'asp__L_c',
                'asn': 'asn__L_c',
                'glu': 'glu__L_c',
                'gln': 'gln__L_c',
                'ser': 'ser__L_c',
                'thr': 'thr__L_c',
                'met': 'met__L_c',
                'cys': 'cys__L_c',
                'lys': 'lys__L_c',
                'arg': 'arg__L_c',
                'his': 'his__L_c',
                'pro': 'pro__L_c',
                'phe': 'phe__L_c',
                'tyr': 'tyr__L_c',
                'trp': 'trp__L_c',
                'X': 'X_c',
                '3-me-glu': '3mglutr_c',  # OK
                'aad': 'L2aadp_c',  # OK
                '4ppro': '4ppro_c',
                'abu': '2abu_c',
                'aeo': '2a89od_c',
                'ala-b': 'ala_B_c',  # OK
                'ala-d': 'ala__D',  # OK
                'allo-thr': 'athr__L_c',  # OK
                'b-ala': 'ala_B_c',
                'beta-ala': 'ala_B_c',
                'bmt': '4b4mt_c',
                'cap': 'cap_c',
                'bht': 'bht_c',
                'dab': '24dab_c',  # OK
                'dhb': 'd23hb_c',  # OK
                'dhpg': '35dhpg_c',
                'dht': 'dht_c',
                'dpg': '35dhpg_c',
                'hiv': '2hiv_c',  # OK
                'hiv-d': 'd2hoiv_c',
                'hmp-d': '2h3mp_c',  # OK
                'horn': 'n5horn_c',  # OK?
                'hpg': '4hpg_c',  # OK ish
                'hyv': '4hylv_c',
                'hyv-d': '2hylv_c',
                'iva': 'iva_c',  # OK
                'lys-b': '36dahx_c',  # OK
                'orn': 'orn_c',  # OK
                'phg': 'phg_c',
                'pip': 'Lpipecol_c',  # OK
                'sal': 'salc_c',  # OK
                'tcl': '555tcl_c',
                'vol': 'valinol',
                'LDAP': '26dap__M_c',  # OK
                'meval': 'meval_c',
                'alle': '2h3mp_c',
                'alaninol': 'alaninol_c',
                'N-(1,1-dimethyl-1-allyl)Trp': 'n11d1at_c',
                'd-lyserg': '99dlys_c',
                'ser-thr': 'thr__L_c',  # OK (is actually serine or threonine)
                'mephe': '99cmphe_c',
                'haorn': 'haorn_c',
                'hasn': 'hasn_c',
                'hforn': 'hforn_c',
                's-nmethoxy-trp': '99snmet_c',
                'alpha-hydroxy-isocaproic-acid': '99ahia_c',
                'MeHOval': '99mehoval',
                '2-oxo-isovaleric-acid': '992oia_c',
                'aoda': '99aoda_c'
                }

t1pks_extenders = {'Malonyl-CoA': 'C00083',
                   'Methylmalonyl-CoA': 'C00683',  # (S)-methylmalonyl-coa. (R) is C01213
                   'Methoxymalonyl-CoA': 'NOKEGG',
                   'Ethylmalonyl-CoA': 'C18026',  # (S)-ethylmalonyl-coa. (R) is C20238
                   'Isobutyryl-CoA': 'C00630',
                   '2-Methylbutyryl-CoA': 'C01033',
                   'trans-1,2-CPDA': 'NOKEGG',
                   'Acetyl-CoA': 'C00024',
                   'Benzoyl-CoA': 'C00512',
                   'Propionyl-CoA': 'C00100',
                   '3-Methylbutyryl-CoA': 'NOKEGG',
                   'CE-Malonyl-CoA': 'NOKEGG',
                   '2-Rhyd-Malonyl-CoA': 'NOKEGG',
                   'CHC-CoA': 'NOKEGG',
                   'inactive': 'NOKEGG',
                   'gly': 'C00037',
                   'ala': 'C00041',
                   'val': 'C00183',
                   'leu': 'C00123',
                   'ile': 'C00407',
                   'asp': 'C00049',
                   'asn': 'C00152',
                   'glu': 'C00025',
                   'gln': 'C00064',
                   'ser': 'C00065',
                   'thr': 'C00188',
                   'met': 'C00073',
                   'cys': 'C00097',
                   'lys': 'C00047',
                   'arg': 'C00062',
                   'his': 'C00135',
                   'pro': 'C00148',
                   'phe': 'C00079',
                   'tyr': 'C00082',
                   'trp': 'C00078',
                   'X': 'CXXXXX'
                   }

std_aa_dic = {
    'gly': 'C00037',
    'ala': 'C00041',
    'val': 'C00183',
    'leu': 'C00123',
    'ile': 'C00407',
    'asp': 'C00049',
    'asn': 'C00152',
    'glu': 'C00025',
    'gln': 'C00064',
    'ser': 'C00065',
    'thr': 'C00188',
    'met': 'C00073',
    'cys': 'C00097',
    'lys': 'C00047',
    'arg': 'C00062',
    'his': 'C00135',
    'pro': 'C00148',
    'phe': 'C00079',
    'tyr': 'C00082',
    'trp': 'C00078',
    'X': 'CXXXXX',
    '3-me-glu': '3mglutr_c',  # OK
    'aad': '2-amino-adipic acid',
    '4ppro': '4-propyl-proline',
    'abu': '2-amino-butyric acid',
    'aeo': '2-amino-8-oxo-9,10-decanoate',
    'ala-b': 'ala_B_c',  # OK
    'ala-d': 'ala__D',  # OK
    'allo-thr': 'athr__L_c',  # OK
    'b-ala': 'ala_B_c',
    'beta-ala': 'ala_B_c',
    'bmt': '# 4-butenyl-4-methyl threonine',
    'cap': 'capreomycidine',
    'bht': 'beta-hydroxy-tyrosine',
    'dab': '24dab_c',  # OK
    'dhb': 'd23hb_c',  # OK
    'dhpg': '3,5-dihydroxy-phenyl-glycine',
    'dht': 'dehydro-threonine/2,3-dehydroaminobutyric acid',
    'dpg': '3,5-dihydroxy-phenyl-glycine (duplicate entry of dhpg)',
    'hiv': '2hiv_c',  # OK
    'hiv-d': 'D-2-hydroxyisovalerate',
    'hmp-d': '2h3mp_c',  # OK
    'horn': 'n5horn_c',  # OK?
    'hpg': '4-hydroxy-phenyl-glycine',
    'hyv': '4-hydroxy-L-valine',
    'hyv-d': '2-hydroxy-valeric acid',
    'iva': '# isovaline',
    'lys-b': '36dahx_c',  # OK
    'orn': 'orn_c',  # OK
    'phg': 'phenyl-glycine',
    'pip': 'Lpipecol_c',  # OK
    'sal': 'salc_c',  # OK
    'tcl': '(4S)-5,5,5-trichloro-leucine',
    'vol': 'valinol',
    'LDAP': '26dap__M_c',  # OK
    'meval': 'Me-Val',
    'alle': '2h3mp_c',
    'alaninol': ' ',
    'N-(1,1-dimethyl-1-allyl)Trp': '',
    'd-lyserg': 'D-lysergic acid',
    'ser-thr': 'thr__L_c',  # OK (is actually serine or threonine)
    'mephe': 'Cmethyl-phenylalanine',
    'haorn': 'L-δ-N-acetyl-δ-N-hydroxyornithine/L-Nδ-hydroxy-Nδ-acylornithine',
    'hasn': 'hydroxyasparagine',
    'hforn': 'L-Nδ-hydroxy-Nδ-formylornithine',
    's-nmethoxy-trp': '',
    'alpha-hydroxy-isocaproic-acid': '',
    'MeHOval': '3-Methyl-2-oxovaleric acid',
    '2-oxo-isovaleric-acid': '',
    'aoda': 'S-2-amino-8-oxodecanoic acid'
}
one_letter_aa = {
    'G': 'gly',
    'P': 'pro',
    'A': 'ala',
    'V': 'val',
    'L': 'leu',
    'I': 'ile',
    'M': 'met',
    'C': 'cys',
    'F': 'phe',
    'Y': 'tyr',
    'W': 'trp',
    'H': 'his',
    'K': 'lys',
    'R': 'arg',
    'Q': 'gln',
    'N': 'asn',
    'E': 'glu',
    'D': 'asp',
    'S': 'ser',
    'T': 'thr'
}
aa_formula = {
    'A': 'C3H7NO2',
    'R': 'C6H14N4O2',
    'N': 'C4H8N2O3',
    'D': 'C4H7NO4',
    'C': 'C3H7NO2S',
    'E': 'C5H9NO4',
    'Q': 'C5H10N2O3',
    'G': 'C2H5NO2',
    'H': 'C6H9N3O2',
    'I': 'C6H13NO2',
    'L': 'C6H13NO2',
    'K': 'C6H14N2O2',
    'M': 'C5H11NO2S',
    'F': 'C9H11NO2',
    'P': 'C5H9NO2',
    'S': 'C3H7NO3',
    'T': 'C4H9NO3',
    'W': 'C11H12N2O2',
    'Y': 'C9H11NO3',
    'V': 'C5H11NO2'
}

acp_domains = ["ACP", "PP-binding", 'ACP_beta']
at_domains = ['Trans-AT_docking', 'PKS_AT']
trans_at_cores = ['transAT-PKS-like', 'transAT-PKS']
nrps_pks_cores = ['T1PKS', 'NRPS', 'PKS-like', 'NRPS-like']
compatible_cores = ['transAT-PKS-like', 'transAT-PKS', 'T1PKS', 'NRPS', 'PKS-like', 'NRPS-like']
RiPPs = ['lanthipeptide', 'thiopeptide', 'lassopeptide']
prefer_types = ['transAT-PKS', 'transAT-PKS-like', 'NRPS', 'T1PKS', 'NRPS-like', 'PKS-like']
dh_domains = ['PKS_DH', 'PKS_DH2']
print_list = ['PKS_KS', 'PKS_DH', 'PKS_DH2', 'PKS_KR', 'MT', 'PKS_ER', 'Condensation', 'Thioesterase']
alternate_starters = ['CAL_domain', 'FkbH', 'GNAT']


def get_gb_list_from_antismash_output(cluster_path):  # yes
    gb_list = []
    for gb_record in SeqIO.parse(open(cluster_path, "r"), "genbank"):
        gb_list.append(gb_record)
    return gb_list

    # Yes this is kinda stupid because i could have just structured this as a dict from the
    # beginning. However, maybe it makes the information extraction process more clear?? :):):)
    # domains, pfam_entries, smCOG, EC_number, rxn_in, rxn_out, core_gene, amino_acid_sequence, predicted_EC, length
    # triple pound signs means that information is stored within a custom class object.


def merge_cores(core_list, high, low):
    '''
    see merge core list
    '''

    for typer in prefer_types:
        for core in core_list:
            if typer == core['type']:
                return [{'start': low, 'end': high,
                         'type': typer,
                         'core_number': core['core_number']}]
            elif core['type'] in RiPPs:
                return [core]

    return [core]


def merge_core_list(core_list):
    '''

    :param core_list:  list of the different core types that have been detected by antiSMASH
                       Have chosen to set boundaries of the core to be equal to the outermost boundaries
                       that antiSMASH predicts. The reason for this is that there are many interleaved clusters
                       bedause they contain both NRPS and PKS modules. this essentially assumes that
                       all core genes are colinear.


    :return:           a single core of type determined by a hierarchy following ->
                       'transAT-PKS', 'transAT-PKS-like', 'NRPS', 'T1PKS', 'NRPS-like', 'PKS-like'
                       transAT PKS can be treated as NRPS and T1PKS, but T1PKS cannot be treated as transAT-PKS
    '''
    low = 1000000000
    high = 0
    for core in core_list:
        if core['start'] < low:
            low = core['start']
        if core['end'] > high:
            high = core['end']

    reduced_list = merge_cores(core_list, high, low)

    return reduced_list


def find_cores_in_cluster(gb_list):
    '''
    In some files, the CDS shows up before the proto core, so tou actually have to find proto_core on beforehand
    (this is in order to be able to know which genes are part of the PKS synthesis)
    This is the reason that we iterate through the gb_list one time before the big next iteration
    :param gb_list:
    :return: A list of cores (cores are overlapping. we are simply interested in the core genes for this project )
    '''
    core_list = []
    backup_core_list = []
    '''
    backup_core_list is there because some BGCs do not have any locations for their clusters. 
    this i think is just for older BGCs that are in mibig, as this does not always happen.
    '''
    core_number = 0
    for gb_record in gb_list:
        for feat in gb_record.features:
            if feat.type == 'source':
                start = feat.location.start
                end = feat.location.end
            if feat.type == 'proto_core':
                core_number += 1

                try:
                    core_list.append(
                        {'start': feat.location.start, 'end': feat.location.end,
                         'type': feat.qualifiers.get('product')[0],
                         'core_number': str(core_number)})
                except AttributeError:

                    backup_core_list.append((
                        {'start': 0, 'end': 0, 'type': feat.qualifiers.get('product')[0],
                         'core_number': str(core_number)}))
                    warnings.warn('cluster has a core with no location')

    try:
        if core_list == []:  # if all cores had no location or no cores exist:
            backup_core_list[0]['start'] = start
            backup_core_list[0]['end'] = end
            return backup_core_list
        else:
            return core_list + backup_core_list
    except IndexError:
        # if the error is due to there being no backup cpre list
        return [{'start': 0, 'end': 0, 'type': 'none',
                 'core_number': '1'}]
    except UnboundLocalError:

        return [{'start': 0, 'end': 0, 'type': 'none',
                 'core_number': '1'}]


def structure_gbk_information(core_list, gb_list):
    CDS = []
    core_domains = {}
    CDS_number = -1
    for gb_record in gb_list:
        for feat in gb_record.features:
            if feat.type == 'CDS':
                if feat.qualifiers.get('gene_functions'):
                    smcog = []
                    gene_function = feat.qualifiers.get('gene_functions')
                    for single in gene_function:
                        if 'SMCOG' in single:
                            smcog.append(single.split('SMCOG')[1].split(':')[0])

                else:
                    smcog = ['No_smCOG']
                cds_is_core = []
                for core in core_list:  # this is wrong, because it adds two instead of 1, every time
                    try:
                        if feat.location.start >= core['start'] and feat.location.end <= core['end']:
                            cds_is_core.append(core)
                    except AttributeError:
                        warnings.warn('Core has no location. Unsure of consequences')
                    # Instantiate CDS-object, consisting of:
                    # list of domains, gene_ontologies, smCOG, EC numbers (possible with many because of many GOs, rxn_in,
                    # rxn_out, core_gene):
                if len(cds_is_core) > 0:  # test to see if the cds is a core gene
                    try:
                        CDS.append(
                            {'smcog': smcog, 'core_gene': cds_is_core, 'domains': [], 'strand': feat.location.strand})
                    except AttributeError:
                        warnings.warn('gene had no location')
                        CDS_number -= 1
                else:
                    try:
                        CDS.append({'smcog': smcog, 'core_gene': False, 'domains': [], 'strand': feat.location.strand})
                    except AttributeError:
                        warnings.warn('gene had no location')
                        CDS_number -= 1
                CDS_number += 1

            else:
                if feat.type == 'aSModule':
                    if feat.qualifiers.get('monomer_pairings'):
                        aids = feat.qualifiers.get('monomer_pairings')
                        for core in core_list:
                            if feat.location.start >= core['start'] and feat.location.end <= core['end']:
                                if core['core_number'] not in core_domains:
                                    core_domains[core['core_number']] = {'type': core['type'], 'modules': [
                                        {'extender_unit': copy.copy(aids)[0].replace('&gt;', '>'),
                                         'start': feat.location.start,
                                         'end': feat.location.end}]}
                                    #  replace('&gt;', '>') -> sometimes there is some parsing error that cant read '>'
                                else:
                                    core_domains[core['core_number']]['modules'] += [
                                        {'extender_unit': copy.copy(aids)[0].replace('&gt;', '>'),
                                         'start': feat.location.start,
                                         'end': feat.location.end}]

                if feat.type == 'CDS_motif':
                    if feat.qualifiers.get('core_sequence'):
                        core_domains[len(core_domains) + 1] = {'type': feat.qualifiers.get('peptide')[0],
                                                               'RiPP': feat.qualifiers.get('core_sequence')[0]}
                if feat.type == 'aSDomain':
                    if feat.qualifiers.get('aSDomain') == ['PKS_AT']:  # AT Domains
                        CDS[CDS_number]['domains'].append(
                            {'type': feat.qualifiers.get('aSDomain')[0],
                             'activity': '',
                             'start': feat.location.start,
                             'end': feat.location.end,
                             'gene': feat.qualifiers.get('locus_tag'),
                             'core_gene': cds_is_core,
                             'minowa': 'mal',
                             'AT_specificity': 'mal'}
                        )
                        for AT_prediction in feat.qualifiers.get('specificity'):
                            if 'consensus: ' in AT_prediction:  # there should only be one instance of this happening
                                CDS[CDS_number]['domains'][-1]['activity'] = AT_prediction.split('consensus: ')[1]
                            elif 'Minowa: ' in AT_prediction:
                                CDS[CDS_number]['domains'][-1]['minowa'] = AT_prediction.split('Minowa: ')[1]
                            elif 'AT_specificity: ' in AT_prediction:
                                CDS[CDS_number]['domains'][-1]['AT_specificity'] = \
                                AT_prediction.split('AT_specificity: ')[1]

                    if feat.qualifiers.get('aSDomain') == ['PKS_KR']:  # KR domains
                        kr_list = []
                        for KR_prediction in feat.qualifiers.get('specificity'):
                            kr_list.append(KR_prediction)
                        for list_item in kr_list:
                            if 'KR activity: ' in list_item:
                                kr_activity = list_item.split('KR activity: ')[1]

                                if kr_activity == 'active':
                                    CDS[CDS_number]['domains'].append(
                                        {'type': feat.qualifiers.get('aSDomain')[0],
                                         'activity': True,
                                         'start': feat.location.start,
                                         'end': feat.location.end,
                                         'gene': feat.qualifiers.get('locus_tag'),
                                         'core_gene': cds_is_core}
                                    )
                                elif kr_activity == 'inactive':
                                    CDS[CDS_number]['domains'].append(
                                        {'type': feat.qualifiers.get('aSDomain')[0],
                                         'activity': False,
                                         'start': feat.location.start,
                                         'end': feat.location.end,
                                         'gene': feat.qualifiers.get('locus_tag'),
                                         'core_gene': cds_is_core}
                                    )
                    elif feat.qualifiers.get('aSDomain') == ['MT']:
                        CDS[CDS_number]['domains'].append(
                            {'type': feat.qualifiers.get('domain_subtype')[0],
                             'activity': True,
                             'start': feat.location.start,
                             'end': feat.location.end,
                             'gene': feat.qualifiers.get('locus_tag'),
                             'core_gene': cds_is_core}
                        )
                    elif feat.qualifiers.get('aSDomain') == ['CAL_domain']:
                        for CAL_prediction in feat.qualifiers.get('specificity'):
                            if 'Minowa: ' in CAL_prediction:
                                CDS[CDS_number]['domains'].append(
                                    {'type': feat.qualifiers.get('aSDomain')[0],
                                     'activity': CAL_prediction.split('Minowa: ')[1],
                                     'start': feat.location.start,
                                     'end': feat.location.end,
                                     'gene': feat.qualifiers.get('locus_tag'),
                                     'core_gene': cds_is_core}
                                )
                    elif feat.qualifiers.get('aSDomain') != ['PKS_AT']:
                        try:
                            CDS[CDS_number]['domains'].append(
                                {'type': feat.qualifiers.get('aSDomain')[0],
                                 'activity': True,
                                 'start': feat.location.start,
                                 'end': feat.location.end,
                                 'gene': feat.qualifiers.get('locus_tag'),
                                 'core_gene': cds_is_core}
                            )
                        except IndexError:
                            # We end up here if the data stored in the gbk-file is weird.
                            # Was necessary in a few cases when the entire mibig database was checked
                            # is probably not relevant when running this code on more recent antismash output
                            print(CDS)
    return {'core_structure': core_domains, 'data': CDS}


def create_json_1(gbk_path):
    gb_list = get_gb_list_from_antismash_output(gbk_path)
    core_list = find_cores_in_cluster(gb_list)

    core_list = merge_core_list(core_list)
    CDS = structure_gbk_information(core_list, gb_list)

    with open(json_path, 'w+') as f:
        json.dump(CDS, f)
    f.close()


def return_key_by_value(eu):
    '''

    :param eu: eu = extender unit. this is found from the gbk file on the format "mal -> ccmal". to get the bigg
    ID of the extender unit, the "mal" part of the string is looked up in one of the lists at the top

    :return: bigg id of the extender unit

    '''

    for long_eu_name in long_to_short:
        if long_to_short[long_eu_name] == eu:
            return long_eu_name
    if eu in std_aa_dic:
        return eu

    return 'pk'


def unfuck_transat_modules(domain_or_module):
    '''

    :param domain_or_module: a list of domains and modules that antismash have found. This is the first of several
    "fixings" of this and only applies to transAS BGCs.

    :return: domain_or_module, except some modules have been broken down if they are ks-
    this might be a truly shit function
    '''
    res = domain_or_module
    for index in range(len(res) - 1):
        if res[index]['element'] == 'module' and res[index + 1]['element'] == 'domain':
            start = index
            if res[start + 1]['info'][
                'type'] in reducing_domains:  # COULD INCLUDE ACP HERE, BUT THOSE CASES MIGHT BE TOO FUCKED
                for domain_loc in range(0, len(domain_or_module[index]['domains']) - 1):
                    # insert domains at end, pop at start
                    if res[index]['domains'][domain_loc]['type'] in acp_domains and \
                            res[index]['domains'][domain_loc + 1]['type'] in at_domains:
                        res = unfuck_transat_modules(res[0:start] + [{'element': 'domain', 'info': x} for x in
                                                                     domain_or_module[index]['domains'][
                                                                     :domain_loc]] + [{'element': 'domain', 'info': x}
                                                                                      for x in domain_or_module[index][
                                                                                                   'domains'][
                                                                                               domain_loc + 1:]] + res[
                                                                                                                   start + 1:])
                return res
    return res


def find_and_replace_DHD_domains(domain_or_module):
    '''
    
    :param domain_or_module: sequence of domains and modules
    :return: find all domain sequences that are KS-DH-ACP and span across two genes. 
    These we put into a non-extending module with extender_unit: 'DHD'
    '''
    res = domain_or_module
    for index in range(len(res)):
        if index + 2 < len(res):
            if res[index]['element'] == res[index + 1]['element'] == res[index + 1]['element'] == 'domain':
                if res[index]['info']['type'] == 'PKS_KS':
                    if res[index + 1]['info']['type'] in dh_domains and res[index + 2]['info']['type'] in acp_domains \
                            and res[index]['info']['gene'] != res[index + 1]['info'][
                        'gene']:  # check that KS and DH is actually on two different genes
                        res = res[0:index] + [{'element': 'module',
                                               'info': {'extender_unit': 'DHD', 'start': res[index]['info']['start'],
                                                        'end': res[index]['info']['end']},
                                               'domains': [res[index]['info'], res[index + 1]['info'],
                                                           res[index + 2]['info']]}] + res[index + 3:]
    return res


def deactivate_reducing_domains_when_kr_is_inactive(domain_or_module):
    '''

    :param domain_or_module:
    :return: res (the same as the input, however this time the er and dh domains have been inactivated if the KR domain
    is inactive. The reason is that antismash only predicts the activity of the KR domains, so by default the DH
    and ER domains are always active. This is impossible if the KR domain is inactive.)
    '''
    res = domain_or_module

    for domod in res:  # domod is a domain or a module
        kr_activity = False
        if domod['element'] == 'module':
            for domain in domod['domains']:
                if domain['type'] == 'PKS_KR':
                    kr_activity = domain['activity']
            for domain in domod['domains']:
                if domain['type'] in dh_er_domains:
                    domain['activity'] = kr_activity

    return res


def activate_kr_domains(domain_or_module):
    '''

    :param domain_or_module:
    :return: for transats, KR domains are rarely inactive, relative to cis-at pks
    therefore the prediction that antismash gives, will in most cases give a worse result than activating all of them
    '''

    res = domain_or_module
    for domod in res:  # domod is a domain or a module
        if domod['element'] == 'module':
            for domain in domod['domains']:
                if domain['type'] == 'PKS_KR':
                    domain['activity'] = True

    return res


def find_and_replace_cut_transat_modules(domain_or_module):
    res = domain_or_module
    module_start = False
    module_end = False
    for index in range(len(res)):
        if domain_or_module[index]['element'] == 'module':
            module_start = False
            module_end = False
        else:
            if domain_or_module[index]['info']['type'] == 'PKS_KS':
                start = index
                module_start = True
            if domain_or_module[index]['info']['type'] in acp_domains and module_start:
                module_end = True
                end = index
            if module_start and module_end:
                res = find_and_replace_cut_transat_modules(res[0:start] + [{'element': 'module',
                                                                            'info': {
                                                                                'extender_unit': 'mal',
                                                                                'start': res[start]['info']['start'],
                                                                                'end': res[end]['info']['end']},
                                                                            'domains': [domain_or_module[a]['info'] for
                                                                                        a in range(start, end)]}] + res[
                                                                                                                    end + 1:])
                return res
    return res


def add_rogue_transat_domains(domain_or_module):
    res = domain_or_module
    # removes duplicate acp-domains - - -
    # This phenoenon is commonly observed for transAT-PKS and from anecdotal evidence they never have any effect on
    # the biosynthesis of the secondary metabolite. The reason they are removed is that this makes it easier to
    # look for modules that are wrongly predicted by antismash

    for ind, domain in enumerate(res):
        if domain['element'] == 'domain':
            if domain['info']['type'] in acp_domains:
                res.pop(ind)

    for index in range(1, len(res) - 1):
        if res[index]['element'] == 'module' and res[index + 1]['element'] == 'domain':

            if res[index + 1]['info']['type'] in reducing_domains and not res[index]['info']['extender_unit'] == 'DHD':
                # case where domain is located right next to module
                not_domain_exists_in_prev_mod = True
                for donk in res[index]['domains']:
                    if res[index + 1]['info']['type'] == donk['type']:
                        # check to see if the rogue domain exists in the module we want to add it to
                        not_domain_exists_in_prev_mod = False
                if not_domain_exists_in_prev_mod:
                    res[index]['domains'].append(res[index + 1]['info'])

            elif res[index + 1]['info']['type'] in reducing_domains and res[index]['info']['extender_unit'] == 'DHD' and \
                    res[index - 1]['element'] == 'module':
                # case where domain is located right next to DHD module, located right next to another module
                not_domain_exists_in_prev_mod = True
                for donk in res[index - 1]['domains']:
                    if res[index + 1]['info']['type'] == donk['type']:
                        # check to see if the rogue domain exists in the module we want to add it to
                        not_domain_exists_in_prev_mod = False

                if not_domain_exists_in_prev_mod:
                    res[index - 1]['domains'].append(res[index + 1]['info'])

    for ind, a in enumerate(res):
        if ind < len(res) - 1:
            if a['element'] == 'module':
                for dom in a['domains']:
                    if dom == res[ind + 1]['info']:
                        res.pop(ind + 1)
                    elif ind < len(res) - 2:
                        if dom == res[ind + 2]['info']:
                            res.pop(ind + 2)

    return res


def add_transat_metabolic_pathway(data, model, core_number, tailoring_reactions):
    """
    :param data: all data that has been structured from gbk-file
    :return: modules max: The maximum detected modules, following rule that it needs a KS-domain and an ACP domain
    Modules mid: All modules that contain KS, AT and ACP domain.
    the modules_min is
    """
    domains = []
    actual_domains = []
    domains_x_modules = []
    reverse_or_positive = 0
    for gene in range(len(data['data'])):
        if data['data'][gene]['core_gene']:
            reverse_or_positive += int(data['data'][gene]['strand'])
    if reverse_or_positive < 0:
        for gene in range(len(data['data'])):
            domains.append([])
            for domain in data['data'][gene]['domains']:
                domains[gene].insert(0, domain)
        domains.reverse()

    else:
        for gene in range(len(data['data'])):
            domains.append([])
            for domain in data['data'][gene]['domains']:
                domains[gene].append(domain)

    for list_of_domains in domains:
        actual_domains += list_of_domains

    for domain in actual_domains:
        domain_is_in_a_module = False
        for module_index, module in enumerate(data['core_structure'][core_number]['modules']):
            if module['start'] <= domain['start'] <= domain['end'] <= module['end']:
                domain_is_in_a_module = True
                if {'element': 'module', 'info': module, 'domains': []} not in domains_x_modules:
                    domains_x_modules.append({'element': 'module', 'info': module, 'domains': []})
        if domain_is_in_a_module == False:
            domains_x_modules.append({'element': 'domain', 'info': domain})

    for domain in actual_domains:
        for module_index, module in enumerate(domains_x_modules):
            if module['info']['start'] <= domain['start'] <= domain['end'] <= module['info']['end'] and module[
                'element'] == 'module':
                domains_x_modules[module_index]['domains'].append(domain)

    # assume that DHD-domains are inactive, however, may still dehydratase the previous domain(?)
    # assue that all others are active.

    domains_x_modules = unfuck_transat_modules(domains_x_modules)
    domains_x_modules = find_and_replace_DHD_domains(domains_x_modules)
    domains_x_modules = find_and_replace_load_modules(domains_x_modules)
    domains_x_modules = find_and_replace_cut_transat_modules(domains_x_modules)
    domains_x_modules = add_rogue_transat_domains(domains_x_modules)

    # domains_x_modules = deactivate_reducing_domains_when_kr_is_inactive(domains_x_modules)
    # the above line seems to be much more untrue than true KR domains are basically always active
    # Thus:
    domains_x_modules = activate_kr_domains(domains_x_modules)

    # The line below is for the case that there is no KR domain present.
    # Because then we dont want dh/er domains to be active
    domains_x_modules = deactivate_reducing_domains_when_kr_is_inactive(domains_x_modules)
    domains_x_modules = search_and_destroy_omt_modules(domains_x_modules)

    # this needs to be last so that we dont add omt domains first as this would delete the previous model.
    # these "rogue" omt domains are usually for tailoring reactions of the pks and is therefore just a bonus
    # It is a weird case, but should make it more correct for the most part, although the methodology behind
    # it is completely wrong.
    model = create_t1_transat_nrps_model(data['core_structure'][core_number], domains_x_modules, model,
                                         tailoring_reactions)

    return model
    # domains x modules is now a dict sorted by placement on BGC.
    # do this to find out if domains are found in this sequence: ks-at-acp


def find_and_replace_load_modules(domain_or_module):
    '''
    :param domain_or_module: dictionary of modules and domains
    :return: domain or module, but add a starter module.

    first check if we have CAL,FkbH or GNAT domains present
    the way this is done is a little wonky, because there aren't really any good ways to do this due to
    how antiSMASH determines the cores.
    we are just interested in the first domains on a gene.
    when we encounter a new gene,
    we want to stop checking if the gene has a starter module.
    we also want to start checking again when we go over to a new gene.

    First we check to see if antismash has located a typical A/AT-ACP/PCP-module:
    if this is the case, we dont want to add a loader module.
    '''
    for module in domain_or_module:
        if module['element'] == 'module':
            if len(module['domains']) == 2:
                if module['domains'][0]['type'] in loader_at_domains and module['domains'][1][
                    'type'] in loader_acp_domains:
                    # (if the first domain is an at/a domain and the second domain is an ACP/PCP-domain:
                    # this solution is a bit hacky, but it should never go through incorrectly based on antiSMASH rules
                    return domain_or_module

    found_loader = False
    prev_gene = ''
    encountered_a_module = False
    for index, dom_mod in enumerate(domain_or_module):  # dom_mod is either a domain or a module
        if dom_mod['element'] == 'domain':  # if this is a domain and not a module
            if not dom_mod['info']['gene'] == prev_gene:
                encountered_a_module = False
                start_of_gene = index  # so we know what domains to remove in case we find a loader domain
                # if we are on a new gene, the search continues. # this is only relevant to do on domains, as modules cannot (should not) by antismash definition# should always be extending, and should in theory not contain cal,fkbh or gnat domains
            if dom_mod['info']['core_gene'] and not encountered_a_module:
                # if this is a core gene. needs to be nested because modules have different structure than domains here
                if dom_mod['info'][
                    'type'] in loader_domains:  # if this is a gene that is a giveaway that this is the loader module:
                    # we have found a loader domain on one of the core genes, before any other modules have been found.
                    found_loader = True
                    loader_type = dom_mod['info']['type']
                    loader_activity = dom_mod['info']['activity']
                    loader_start = start_of_gene
                    loader_gene = dom_mod['info']['gene']
            prev_gene = dom_mod['info']['gene']
            if found_loader and dom_mod['info']['type'] in loader_acp_domains:
                return domain_or_module[0:start_of_gene] + [{'element': 'module', 'info': {
                    'extender_unit': loader_type,
                    'loader_activity': loader_activity,
                    'start': domain_or_module[start_of_gene]['info']['start'],
                    'end': domain_or_module[index]['info']['end']}, 'domains': [x['info'] for x in domain_or_module[
                                                                                                   loader_start:index]]}] + domain_or_module[
                                                                                                                            index:]

        elif dom_mod['element'] == 'module':
            encountered_a_module = True  # to check that we have not passed any real modules
            if found_loader and dom_mod['domains'][0]['gene'] == prev_gene:
                return domain_or_module[0:start_of_gene] + [{'element': 'module', 'info': {
                    'extender_unit': loader_type,
                    'loader_activity': loader_activity,
                    'start': domain_or_module[start_of_gene]['info']['start'],
                    'end': domain_or_module[index]['info']['end']}, 'domains': [x['info'] for x in domain_or_module[
                                                                                                   loader_start:index - 1]]}] + domain_or_module[
                                                                                                                                index:]

            if not dom_mod['domains'][0]['gene'] == prev_gene:
                start_of_gene = index  # so we know what domains to remove in case we find a loader domain
            prev_gene = dom_mod['domains'][0]['gene']
            '''
            the [0] part of the above line is just because modules can in theory span across several genes
            '''

    # Then we check if there are any loader modules that can be found within the core genes that are not inside a
    # module as defined by antiSMASH. These are the [A/AT]-[ACP/PCP/PP] modules.
    prev_domain = ''
    start_of_loader = 0
    prev_act = ''
    for index, dom_mod in enumerate(domain_or_module):  # dom_mod is either a domain or a module
        if dom_mod['element'] == 'domain':  # if this is a domain and not a module (implicitly not a module)
            if not dom_mod['info']['gene'] == prev_gene:  # If the domain is the first domain on a gene.
                started_new_gene = True

            # starter_module
            if dom_mod['info']['core_gene'] and not encountered_a_module:
                # if this is a core gene. needs to be nested because modules have different structure than domains
                if prev_domain in loader_at_domains and dom_mod['info']['type'] in loader_acp_domains and prev_gene == \
                        dom_mod['info']['gene']:
                    # if this is a gene that contains the [A/AT]-[ACP/PCP/PP] module outside of a regular
                    # extending module:
                    return domain_or_module[0:start_of_loader] + [{'element': 'module', 'info': {
                        'extender_unit': 'starter_unit',
                        'activity': prev_act,
                        'start': domain_or_module[start_of_loader]['info']['start'],
                        'end': domain_or_module[index]['info']['end']},
                                                                   'domains': [domain_or_module[index - 1]['info'],
                                                                               domain_or_module[index][
                                                                                   'info']]}] + domain_or_module[index:]
                start_of_loader = index
            prev_domain = dom_mod['info']['type']
            prev_act = dom_mod['info']['activity']
            prev_gene = dom_mod['info']['gene']


        elif dom_mod['element'] == 'module':
            encountered_a_module = True  # to check that we have not passed any real modules
            prev_gene = dom_mod['domains'][0]['gene']

    # in case we dont find any loader module, we assume its existence anyways, based on the finding that
    # this is usually the case. We first try to find a free-standing A or AT-domain in order to predict wether the
    # Starter is an amino acid or a Acyl-unit.
    domain_or_module.insert(0, {'element': 'module', 'info': {
        'extender_unit': 'custom_starter', 'start': 0,
        'end': 0}, 'domains': []})

    return domain_or_module


def search_and_destroy_omt_modules(domain_or_module):
    # modules with omt domains are often inactive.
    # so we change the extender unit to oMT so we can catch this later
    res = copy.copy(domain_or_module)
    for index, module in enumerate(domain_or_module):
        if module['element'] == 'module':
            for domain in module['domains']:
                if domain['type'] == 'oMT':
                    res[index]['info']['extender_unit'] = 'oMT'
    return res


def add_t1pks_metabolic_pathway(data, model, core_number, tailoring_reactions):
    """
        :param data: all data that has been structured from gbk-file
        :return: modules max: The maximum detected modules, following rule that it needs a KS-domain and an ACP domain
        Modules mid: All modules that contain KS, AT and ACP domain.
        the modules_min is
        """
    domains = []
    actual_domains = []
    domains_x_modules = []
    reverse_or_positive = 0
    for gene in range(len(data['data'])):
        if data['data'][gene]['core_gene']:
            reverse_or_positive += int(data['data'][gene]['strand'])

    if reverse_or_positive < 0:
        for gene in range(len(data['data'])):
            domains.append([])
            for domain in data['data'][gene]['domains']:
                domains[gene].insert(0, domain)
        domains.reverse()
    else:
        for gene in range(len(data['data'])):
            domains.append([])
            for domain in data['data'][gene]['domains']:
                domains[gene].append(domain)

    for list_of_domains in domains:
        actual_domains += list_of_domains

    for domain in actual_domains:
        domain_is_in_a_module = False
        for module_index, module in enumerate(data['core_structure'][core_number]['modules']):
            if module['start'] <= domain['start'] <= domain['end'] <= module['end']:
                domain_is_in_a_module = True
                if {'element': 'module', 'info': module, 'domains': []} not in domains_x_modules:
                    domains_x_modules.append({'element': 'module', 'info': module, 'domains': []})
        if domain_is_in_a_module == False:
            domains_x_modules.append({'element': 'domain', 'info': domain})

    for domain in actual_domains:
        for module_index, module in enumerate(domains_x_modules):
            if module['info']['start'] <= domain['start'] <= domain['end'] <= module['info']['end'] and module[
                'element'] == 'module':
                domains_x_modules[module_index]['domains'].append(domain)

    # assume that DHD-domains are inactive, however, may still dehydratase the previous domain(?)
    # assue that all others are active.

    domains_x_modules = find_and_replace_load_modules(domains_x_modules)
    domains_x_modules = deactivate_reducing_domains_when_kr_is_inactive(domains_x_modules)
    model = create_t1_transat_nrps_model(data['core_structure'][core_number], domains_x_modules, model,
                                         tailoring_reactions)


    return model
    # domains x modules is now a dict sorted by placement on BGC.
    # do this to find out if domains are found in this sequence: ks-at-acp


def add_nrps_metabolic_pathway(data, model, core_number, tailoring_reactions):
    """
        :param data: all data that has been structured from gbk-file
        :return: modules max: The maximum detected modules, following rule that it needs a KS-domain and an ACP domain
        Modules mid: All modules that contain KS, AT and ACP domain.
        the modules_min is
        """
    domains = []
    actual_domains = []
    domains_x_modules = []
    reverse_or_positive = 0
    for gene in range(len(data['data'])):
        if data['data'][gene]['core_gene']:
            reverse_or_positive += int(data['data'][gene]['strand'])
    if reverse_or_positive < 0:
        for gene in range(len(data['data'])):
            domains.append([])
            for domain in data['data'][gene]['domains']:
                domains[gene].insert(0, domain)
        domains.reverse()
    else:
        for gene in range(len(data['data'])):
            domains.append([])
            for domain in data['data'][gene]['domains']:
                domains[gene].append(domain)

    for list_of_domains in domains:
        actual_domains += list_of_domains

    for domain in actual_domains:
        domain_is_in_a_module = False
        for module_index, module in enumerate(data['core_structure'][core_number]['modules']):
            if module['start'] <= domain['start'] <= domain['end'] <= module['end']:
                domain_is_in_a_module = True
                if {'element': 'module', 'info': module, 'domains': []} not in domains_x_modules:
                    domains_x_modules.append({'element': 'module', 'info': module, 'domains': []})
        if domain_is_in_a_module == False:
            domains_x_modules.append({'element': 'domain', 'info': domain})

    for domain in actual_domains:
        for module_index, module in enumerate(domains_x_modules):
            if module['info']['start'] <= domain['start'] <= domain['end'] <= module['info']['end'] and module[
                'element'] == 'module':
                domains_x_modules[module_index]['domains'].append(domain)

    domains_x_modules = find_and_replace_load_modules(domains_x_modules)
    model = create_t1_transat_nrps_model(data['core_structure'][core_number], domains_x_modules, model,
                                         tailoring_reactions)


    return model


def force_X_nrps_module_flux(substrate):
    '''
    :param substrate: an amino acid that is either not specifically determined (NRPS predictor couldnt find
    specific amino acid)
    :return: 20 reactions that convert each amino acid into the unknown amino acid
    Is actually not that terrible of an assumption (albeit with respect to linear program solving which is not that
    good in the case of secondary metabolites, but i digress)
    consider the biosynthesis of HPG which can be the literal conversion of tyrosine to hpg, using only o2 as a
    cofactor.
    '''

    metabolite = cobra.Metabolite(substrate + '_c', formula='na', name='Unknown amino acid', compartment='c')
    reaction_list_aa = []
    for amino_acid in acid_to_bigg:
        reaction = cobra.Reaction(amino_acid + '_to_' + substrate)
        reaction.name = 'Convert known to unknown_amino_acid'
        reaction.lower_bound = 0.  # This is the default
        reaction.upper_bound = 1000.
        reaction.add_metabolites({snorre_model.metabolites.get_by_id(acid_to_bigg[amino_acid]): -1,
                                  metabolite: 1})
        reaction_list_aa.append(reaction)
    return reaction_list_aa


def create_t1_transat_nrps_model(core_structure, domains_x_modules, model, tailoring_reactions):
    # DUMMY-REACTION som er startmetabolitt ettersom det alltid kreves at det finnes en tidligere undermetabolitt når
    # vi lager modellen som inneholder alle reaksjonene.
    reaction_list = []  # list of reactions that take place
    if tailoring_reactions['methoxymalonyl']:
        reaction = cobra.Reaction('Mxmal-ACP')
        reaction.name = 'Methoxymalonyl-ACP synthesis'
        reaction.lower_bound = 0.  # This is the default
        reaction.upper_bound = 1000.
        reaction.add_metabolites(cofactor_reactions_dict['mxmal'])
        reaction_list.append(reaction)

    in_rx = cobra.Reaction('IN_' + core_structure['type'])
    in_rx.name = core_structure['type'] + '_reaction'
    in_rx.lower_bound = 0.  # This is the default
    in_rx.upper_bound = 1000.
    in_rx.add_metabolites({cobra.Metabolite(core_structure['type'] + '_0',
                                            formula='dummy',
                                            name=core_structure['type'],
                                            compartment='d'): 1.0})

    domain_counter = 1

    for module in domains_x_modules:

        if module['element'] == 'module' and not module['info']['extender_unit'] == 'DHD':
            # Little hack to find out if it is a nrps or pks module. PKS_KS and Condensation
            # domains are obligatory for their respective type of module
            # This is necessary because sometimes we find metabolites that are not in the snorre_model
            module_type = ''
            if module['info']['extender_unit'] == 'custom_starter':
                module_type = 'PKS'

            for domain in module['domains']:

                if domain['type'] == 'PKS_KS' or domain['type'] == 'PKS_AT':
                    module_type = 'PKS'
                elif domain['type'] == 'Condensation' or domain['type'] == 'AMP-binding':
                    module_type = 'NRPS'

            # This needs to be first because we want to add extender unit as first step of a module, but for transAT
            # modules, the at domain is not always present.

            if module['info']['extender_unit'] not in exclude_modules:
                # essentially: if the module is a regular extender module
                try:
                    # This is just a way that we can separate extender units that exist in the GEM
                    # from those that do not exist in the GEM. if the TRY fails, it means that
                    # the extender unit does not exist in the GEM
                    if not module_type == 'starter_module':
                        extender_unit = return_key_by_value(module['info']['extender_unit'].split(' -> ')[0])
                    else:
                        extender_unit = return_key_by_value(module['info']['activity'].split(' -> ')[0])
                    # happens if extender unit is not malcoa, mmcoa, mxcoa or emcoa,
                    # (or minowa and at_specificity disagrees):
                    if extender_unit == 'pk':

                        # if this is the empty loader we created:
                        if module['info']['extender_unit'] == 'custom_starter':
                            extender_unit = 'Malonyl-CoA'

                        # if this is a alternate loader, i.e gnat, fkbh, cal loader:
                        elif module['info']['extender_unit'] in alternate_starters:
                            if module['info']['extender_unit'] == 'CAL_domain':
                                module_type = module['info']['loader_activity']
                            elif module['info']['extender_unit'] == 'FkbH':
                                module_type = 'FkbH'
                            elif module['info']['extender_unit'] == 'GNAT':
                                module_type = 'GNAT'
                            elif module['info']['extender_unit'] == 'starter_unit':
                                module_type = 'starter_module'

                        # then  we check if this should be mxmal:
                        # if minowa and at_prediction had no consensus, we set extender unit to mxmal if there is
                        # machinery that synthesises it. if not, we assume the extender is malcoa
                        for domains in module['domains']:
                            if domains['type'] == 'PKS_AT':
                                if tailoring_reactions['methoxymalonyl'] and domains['minowa'] == 'mxmal' or domains[
                                    'AT_specificity'] == 'mxmal':
                                    extender_unit = 'Methoxymalonyl-CoA'
                                else:
                                    # we have to chose between minowa and AT_specificity, so we chose minowa:
                                    # this also has the consequence that if minowa is emal: the next step may convert
                                    # the extender into methoxymalonyl-CoA
                                    extender_unit = return_key_by_value(domains['minowa'])

                    if tailoring_reactions['methoxymalonyl'] and extender_unit == 'Ethylmalonyl-CoA':
                        extender_unit = 'Methoxymalonyl-CoA'
                    prevmet = cobra.Metabolite(core_structure['type'] + '_' + str(domain_counter - 1),
                                               formula='unknown',
                                               name=core_structure['type'],
                                               compartment='d')
                    postmet = cobra.Metabolite(core_structure['type'] + '_' + str(domain_counter),
                                               formula='unknown',
                                               name=core_structure['type'],
                                               compartment='d')

                    reaction = cobra.Reaction(core_structure['type'] + '_' + str(domain_counter))
                    reaction.name = core_structure['type'] + '_reaction_' + str(domain_counter)
                    reaction.lower_bound = 0.  # This is the default
                    reaction.upper_bound = 1000.
                    if module_type == 'NRPS':
                        reaction.add_metabolites(
                            {snorre_model.metabolites.get_by_id(long_to_bigg[extender_unit]): -1,
                             snorre_model.metabolites.get_by_id('atp_c'): -1,
                             snorre_model.metabolites.get_by_id('amp_c'): 1,
                             snorre_model.metabolites.get_by_id('ppi_c'): 1,
                             prevmet: -1,
                             postmet: 1})
                        reaction_list.append(reaction)
                        domain_counter += 1
                    elif module_type == 'PKS':
                        if extender_unit == 'Methoxymalonyl-CoA' or extender_unit == 'Ethylmalonyl-CoA':
                            reaction.add_metabolites(
                                {cobra.Metabolite(long_to_bigg[extender_unit], formula='unknown', name=extender_unit,
                                                  compartment='c'): -1,
                                 snorre_model.metabolites.get_by_id('coa_c'): 1,
                                 snorre_model.metabolites.get_by_id('co2_c'): 1,
                                 prevmet: -1,
                                 postmet: 1})
                        else:
                            reaction.add_metabolites(
                                {snorre_model.metabolites.get_by_id(long_to_bigg[extender_unit]): -1,
                                 snorre_model.metabolites.get_by_id('coa_c'): 1,
                                 snorre_model.metabolites.get_by_id('co2_c'): 1,
                                 prevmet: -1,
                                 postmet: 1})
                        reaction_list.append(reaction)

                        domain_counter += 1


                    elif module_type == 'FkbH':
                        reaction.add_metabolites(
                            {prevmet: -1,
                             postmet: 1})
                        reaction.add_metabolites(cofactor_reactions_dict['fkbh'])
                        reaction_list.append(reaction)
                        domain_counter += 1
                    elif module_type == 'GNAT':
                        reaction.add_metabolites(
                            {prevmet: -1,
                             postmet: 1})
                        reaction.add_metabolites(cofactor_reactions_dict['gnat'])
                        reaction_list.append(reaction)
                        domain_counter += 1
                    elif module_type == 'fatty_acid':
                        reaction.add_metabolites(
                            {prevmet: -1,
                             postmet: 1})
                        reaction.add_metabolites(cofactor_reactions_dict['fatty_acid'])
                        reaction_list.append(reaction)
                        domain_counter += 1
                    elif module_type == 'AHBA':
                        reaction.add_metabolites(
                            {prevmet: -1,
                             postmet: 1})
                        reaction.add_metabolites(cofactor_reactions_dict['ahba'])
                        reaction_list.append(reaction)
                        domain_counter += 1
                    elif module_type == 'shikimic_acid':
                        reaction.add_metabolites(
                            {prevmet: -1,
                             postmet: 1})
                        reaction.add_metabolites(cofactor_reactions_dict['shikimic_acid'])
                        reaction_list.append(reaction)
                        domain_counter += 1
                    elif module_type == 'Acetyl-CoA':
                        reaction.add_metabolites(
                            {prevmet: -1,
                             postmet: 1})
                        reaction.add_metabolites(cofactor_reactions_dict['acetyl'])
                        reaction_list.append(reaction)
                        domain_counter += 1
                    elif module_type == 'NH2':
                        reaction.add_metabolites(
                            {prevmet: -1,
                             postmet: 1})
                        reaction.add_metabolites(cofactor_reactions_dict['NH2'])
                        reaction_list.append(reaction)
                        domain_counter += 1
                    elif module_type == 'starter_module':
                        raise IndexError
                        reaction.add_metabolites(
                            {prevmet: -1,
                             postmet: 1})
                        reaction.add_metabolites(cofactor_reactions_dict['NH2'])
                        reaction_list.append(reaction)
                        domain_counter += 1
                    else:
                        # raise arbitrary error because something has gone terribly wrong
                        raise IndexError
                except KeyError:  # the substrate is not in the model you want to insert into
                    
                    warnings.warn('Predicted substrate does not exist in the model')
                    extender_unit = return_key_by_value(module['info']['extender_unit'].split(' -> ')[0])

                    if module_type == 'PKS':
                        extender_unit = 'Malonyl-CoA'
                        # the only cases we end up here when the module type is PKS, is when a
                        # starter unit is predicted as the extender unit. This Cannot be true for starters that antiSMASH
                        # predicts, so we set this extender to malonyl-CoA

                    prevmet = cobra.Metabolite(core_structure['type'] + '_' + str(domain_counter - 1),
                                               formula='unknown',
                                               name=core_structure['type'],
                                               compartment='c')
                    postmet = cobra.Metabolite(core_structure['type'] + '_' + str(domain_counter),
                                               formula='unknown',
                                               name=core_structure['type'],
                                               compartment='c')
                    reaction = cobra.Reaction(
                        core_structure['type'] + '_' + str(domain_counter))  # this line is different
                    reaction.name = core_structure['type'] + '_reaction_' + str(domain_counter)
                    reaction.lower_bound = 0.  # This is the default
                    reaction.upper_bound = 1000.
                    if module_type == 'NRPS':
                        if extender_unit == 'hpg':
                            # the 4 next reactions enable the synthesis of hpg from common metabolites
                            hpg_reaction1 = cobra.Reaction('hpg_synthesis_1')
                            hpg_reaction1.name = 'synthesis of hpg'
                            hpg_reaction1.lower_bound = 0.  # This is the default
                            hpg_reaction1.upper_bound = 1000.
                            hpg_reaction1.add_metabolites(cofactor_reactions_dict['hpg_1'])

                            hpg_reaction2 = cobra.Reaction('hpg_synthesis_2')
                            hpg_reaction2.name = 'synthesis of hpg'
                            hpg_reaction2.lower_bound = 0.  # This is the default
                            hpg_reaction2.upper_bound = 1000.
                            hpg_reaction2.add_metabolites(cofactor_reactions_dict['hpg_2'])

                            hpg_reaction3 = cobra.Reaction('hpg_synthesis_3')
                            hpg_reaction3.name = 'synthesis of hpg'
                            hpg_reaction3.lower_bound = 0.  # This is the default
                            hpg_reaction3.upper_bound = 1000.
                            hpg_reaction3.add_metabolites(cofactor_reactions_dict['hpg_3'])

                            hpg_reaction4 = cobra.Reaction('hpg_synthesis_4')
                            hpg_reaction4.name = 'synthesis of hpg'
                            hpg_reaction4.lower_bound = 0.  # This is the default
                            hpg_reaction4.upper_bound = 1000.
                            hpg_reaction4.add_metabolites(cofactor_reactions_dict['hpg_4'])

                            model.add_reactions([hpg_reaction1, hpg_reaction2, hpg_reaction3, hpg_reaction4])

                        elif extender_unit == 'dpg' or extender_unit == 'dhpg':
                            # duplicate entries
                            # dpg and dhpg are the same substrates
                            dpg_reaction1 = cobra.Reaction('hpg_synthesis')
                            dpg_reaction1.name = 'synthesis of dhpg'
                            dpg_reaction1.lower_bound = 0.  # This is the default
                            dpg_reaction1.upper_bound = 1000.
                            dpg_reaction1.add_metabolites(cofactor_reactions_dict['dpg'])
                            model.add_reactions([dpg_reaction1])
                        elif extender_unit == 'pip':
                            # duplicate entries
                            # dpg and dhpg are the same substrates
                            pip_reaction1 = cobra.Reaction('pipecolic_acid_synthesis')
                            pip_reaction1.name = 'synthesis of pipecolic acid'
                            pip_reaction1.lower_bound = 0.  # This is the default
                            pip_reaction1.upper_bound = 1000.
                            pip_reaction1.add_metabolites(cofactor_reactions_dict['pip'])
                            model.add_reactions([pip_reaction1])
                        else:
                            # we dont know the specificity of the A-domain, so we let any amino acid take the role
                            # as this amino acid. If you want this to be true for ALL amino acids (i.e. haorn, pip etc.)
                            # change the elif above to a simple "else"
                            model.add_reactions(force_X_nrps_module_flux(extender_unit))
                        reaction.add_metabolites({
                            cobra.Metabolite(long_to_bigg[extender_unit], formula='unknown', name=extender_unit,
                                             compartment='c'): -1,
                            snorre_model.metabolites.get_by_id('atp_c'): -1,
                            snorre_model.metabolites.get_by_id('amp_c'): 1,
                            snorre_model.metabolites.get_by_id('ppi_c'): 1,
                            prevmet: -1,
                            postmet: 1})
                        reaction_list.append(reaction)
                        warnings.warn('Predicted substrate does not exist in the model')
                        domain_counter += 1
                    elif module_type == 'PKS':
                        reaction.add_metabolites(
                            {cobra.Metabolite(long_to_bigg[extender_unit], formula='unknown', name=extender_unit,
                                              compartment='c'): -1,
                             snorre_model.metabolites.get_by_id('coa_c'): 1,
                             snorre_model.metabolites.get_by_id('co2_c'): 1,
                             prevmet: -1,
                             postmet: 1})
                        reaction_list.append(reaction)
                        domain_counter += 1
                    else:
                        raise IndexError

            for domains in module['domains']:
                if domains['activity']:

                    if domains['type'] in general_domain_dict:
                        reaction = cobra.Reaction(core_structure['type'] + '_' + str(domain_counter))
                        reaction.name = core_structure['type'] + '_reaction_' + str(domain_counter)
                        reaction.lower_bound = 0.  # This is the default
                        reaction.upper_bound = 1000.
                        reaction.add_metabolites(cofactor_reactions_dict[general_domain_dict[domains['type']]])
                        prevmet = cobra.Metabolite(core_structure['type'] + '_' + str(domain_counter - 1),
                                                   formula='unknown',
                                                   name=core_structure['type'],
                                                   compartment='d')
                        postmet = cobra.Metabolite(core_structure['type'] + '_' + str(domain_counter),
                                                   formula='unknown',
                                                   name=core_structure['type'],
                                                   compartment='d')
                        reaction.add_metabolites({prevmet: -1, postmet: 1})
                        reaction_list.append(reaction)
                        domain_counter += 1

    '''
    this part below is to add tailoring reactions:
    tailoring reactions is a dict of type:

    {glycosyltransferase: 3, ALA: 0, glycerol: 1}
    '''

    for tailoring_reaction in tailoring_reactions:
        if tailoring_reaction != 'methoxymalonyl':
            if tailoring_reactions[tailoring_reaction]:
                for repetition in range(tailoring_reactions[tailoring_reaction]):
                    reaction = cobra.Reaction(core_structure['type'] + '_' + str(domain_counter))
                    reaction.add_metabolites(tailoring_reactions_dict[tailoring_reaction])
                    prevmet = cobra.Metabolite(core_structure['type'] + '_' + str(domain_counter - 1),
                                               formula='unknown',
                                               name=core_structure['type'],
                                               compartment='d')
                    postmet = cobra.Metabolite(core_structure['type'] + '_' + str(domain_counter),
                                               formula='unknown',
                                               name=core_structure['type'],
                                               compartment='d')
                    reaction.add_metabolites({prevmet: -1.0, postmet: 1.0})
                    reaction_list.append(reaction)
                    domain_counter += 1

    '''
    The reaction below is here in order to have a reaction that converts e.g. T1PKS_54 to a general metabolite
    we do this so that we can add tailoring reactions. 
    '''
    ex_rx = cobra.Reaction('EX_secondary_metabolite')
    ex_rx.name = core_structure['type'] + '_reaction'
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.

    ex_rx.add_metabolites({cobra.Metabolite(core_structure['type'] + '_' + str(domain_counter - 1),
                                            formula='unknown',
                                            name=core_structure['type'],
                                            compartment='d'): - 1.0})

    model.add_reactions(reaction_list + [ex_rx] + [in_rx])
    lump_metabolites = {}
    for rx in reaction_list:
        for reactant in rx.metabolites:
            if not reactant in lump_metabolites:
                lump_metabolites[reactant] = rx.metabolites[reactant]
            else:
                lump_metabolites[reactant] += rx.metabolites[reactant]

    lump_model = cobra.Model('lump_model')
    lump_rx = cobra.Reaction('LUMP_' + core_structure['type'])
    in_rx.name = core_structure['type'] + '_reaction'
    in_rx.lower_bound = 0.  # This is the default
    in_rx.upper_bound = 1000.
    lump_rx.add_metabolites(lump_metabolites)
    lump_model.add_reactions([lump_rx])
    return model, lump_model


def add_ripp_metabolic_pathway(core_structure, model):
    reaction = cobra.Reaction(core_structure['type'])
    reaction.name = core_structure['type'] + '_reaction'
    reaction.lower_bound = 0.  # This is the default
    reaction.upper_bound = 1000.
    reaction_metabolites = {cobra.Metabolite(core_structure['type'],
                                             formula='unknown_but_can_find',
                                             name=core_structure['type'],
                                             compartment='d'): 1.0}
    aa_metabolites = {}
    for letter in core_structure['RiPP']:
        if letter in aa_metabolites:
            aa_metabolites[letter] -= 1
        else:
            aa_metabolites[letter] = -1

    for metabolite in aa_metabolites:  # color red and metabolites_snorre

        try:
            met_bigg = long_to_bigg[one_letter_aa[metabolite]]
            metabolite_added = snorre_model.metabolites.get_by_id(met_bigg)
            reaction_metabolites[metabolite_added] = aa_metabolites[metabolite]
        except KeyError:
            continue

    reaction.add_metabolites(reaction_metabolites)
    model.add_reactions([reaction])

    return model


def find_tailoring_reactions_from_smcogs(data):
    '''

    :param data: the data collected from the gbk file aquired from antiSMASH - here we are interested in the smCOGs
    of genes. genes inside the core genes are excluded, meaning that no false positives from fkbh loaders could lead to
    incorrect assumptions of a methoxymalonyl-coa pathway existing.
    :return: dict containing information on if the specific smCOG gene is found in the cluster.
    This is processed by function "add_tailoring_smcogs" (tailoring reactions) and process_extender_unit_smcogs
    (essentially just to find methoxymalonyl-coa synthesis genes)
    '''

    tailoring_genes_based_on_smCOG_definition_dict = {}

    tailoring_genes_based_on_smCOG_definition_dict[
        1256] = False  # FkbH_like protein with smcog number 1256 - for methoxymalonyl-coa
    tailoring_genes_based_on_smCOG_definition_dict[
        1095] = False  # 3-Hydroxybutyryl-CoA dehydrogenase like protein with smcog number 1095 - for methoxymalonyl-coa
    tailoring_genes_based_on_smCOG_definition_dict[
        1062] = False  # glycosyltransferase with smcog 1062 - For glycosyl groups
    tailoring_genes_based_on_smCOG_definition_dict[
        1084] = False  # for tailoring reaction that adds glycerol to PK through 1,3-biphosphoglycerate
    tailoring_genes_based_on_smCOG_definition_dict[
        1002] = False  # AMP-dependent synthase and ligase with smcog 1002 - 2-Amino-3-hydroxycyclopent-2-enone
    tailoring_genes_based_on_smCOG_definition_dict[
        1109] = False  # 8-amino-7-oxononanoate synthase with smcog 1109- For 2-Amino-3-hydroxycyclopent-2-enone
    for gene in range(len(data['data'])):
        for smcog in data['data'][gene]['smcog']:
            if smcog == '1256':
                tailoring_genes_based_on_smCOG_definition_dict[1256] = True
            elif smcog == '1095':
                tailoring_genes_based_on_smCOG_definition_dict[1095] = True
            elif smcog == '1062':
                if not tailoring_genes_based_on_smCOG_definition_dict[1062]:
                    tailoring_genes_based_on_smCOG_definition_dict[1062] = 1
                else:
                    tailoring_genes_based_on_smCOG_definition_dict[1062] += 1
            elif smcog == '1084':
                tailoring_genes_based_on_smCOG_definition_dict[1084] = True
            elif smcog == '1109':
                tailoring_genes_based_on_smCOG_definition_dict[1109] = True
                # the code beneath is to check genes flanking the 1109-gene if it is an amp-synthase and ligase
                try:
                    for smcog in data['data'][gene - 1]['smcog']:
                        if smcog == '1002':
                            tailoring_genes_based_on_smCOG_definition_dict[1002] = True
                except Exception:
                    continue
                try:
                    for smcog in data['data'][gene + 1]['smcog']:
                        if smcog == '1002':
                            tailoring_genes_based_on_smCOG_definition_dict[1002] = True
                except Exception:
                    continue

    return tailoring_genes_based_on_smCOG_definition_dict


def add_tailoring_smcogs(smcog_dict):  # connect these reactions to the metabolite 'secondary_metabolite'
    res = {'ALA': 0, 'glycosyltransferase': 0, 'glycerol': 0, 'methoxymalonyl': 0}
    if smcog_dict[1002] and smcog_dict[1109]:
        res['ALA'] = 1
    if smcog_dict[1256] and smcog_dict[1084]:
        # add pathway for adding 1,3-biphosphoglycerate
        res['glycerol'] = 1
    if smcog_dict[1062]:
        for a in range(smcog_dict[1062]):
            res['glycosyltransferase'] += 1
    if smcog_dict[1256] and smcog_dict[1095]:
        res['methoxymalonyl'] = 1

        # add number of glycosides equal to smcog_dict[1062]
    return res


def add_cores_to_model(data, model_output_path):
    model = cobra.Model('BGC_name')

    smcog_dict = find_tailoring_reactions_from_smcogs(data)

    tailoring_reactions_dict = add_tailoring_smcogs(smcog_dict)
    '''
    First we check to see if methoxymalonyl-coa is synthesized by the BGC:
    '''

    for core_number in data['core_structure']:
        if data['core_structure'][core_number]['type'] in RiPPs:
            model = add_ripp_metabolic_pathway(data['core_structure'][core_number], model)

        elif data['core_structure'][core_number]['type'] == 'transAT-PKS':
            model, lump_model = add_transat_metabolic_pathway(data, model, core_number, tailoring_reactions_dict)

        elif data['core_structure'][core_number]['type'] == 'T1PKS':
            model, lump_model = add_t1pks_metabolic_pathway(data, model, core_number, tailoring_reactions_dict)

        elif data['core_structure'][core_number]['type'] == 'NRPS':
            model, lump_model = add_nrps_metabolic_pathway(data, model, core_number, tailoring_reactions_dict)
    '''
    Then finally, we save the model to a file
    '''
    cobra.io.save_json_model(model, model_output_path)


if __name__ == '__main__':
    for filename in os.listdir(biggbk):
        if filename.endswith('.gbk'):
            json_path = json_folder + filename.split('.')[0] + '.json'
            create_json_1(biggbk + filename)
            with open(json_path, 'r') as json_file:
                data_json_1 = json.load(json_file)
            json_file.close()
            add_cores_to_model(data_json_1, output_gbk + filename[:-4] + ".json")
