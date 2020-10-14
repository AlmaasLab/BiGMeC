#!/usr/bin/env python
# coding: utf-8
"""
Authors:
  - Snorre Sulheim, snorres.sulheim@sintef.no
  - Fredrik Fossheim

Date: 17.09.2020
Lisence: CC-BY-4.0

This file contains dictionaries that are used to look up metabolite and reactions names, IDs etc.

"""
import cobra
model_fn = '../Models/BiGG_universal_model.json'
ref_model = cobra.io.load_json_model(model_fn)

cofactor_metabolites_dict = {
    'nadph': ref_model.metabolites.get_by_id('nadph_c'),
    'nadp': ref_model.metabolites.get_by_id('nadp_c'),
    'h2o': ref_model.metabolites.get_by_id('h2o_c'),
    'coa': ref_model.metabolites.get_by_id('coa_c'),
    'co2': ref_model.metabolites.get_by_id('co2_c'),
    'h+': ref_model.metabolites.get_by_id('h_c'),
    'sam': ref_model.metabolites.get_by_id('amet_c'),
    'sah': ref_model.metabolites.get_by_id('ahcys_c'),
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
    'mxmal': {ref_model.metabolites.get_by_id('coa_c'): -1.0,
              ref_model.metabolites.get_by_id('13dpg_c'): -1.0,
              ref_model.metabolites.get_by_id('pi_c'): 2.0,
              ref_model.metabolites.get_by_id('nadp_c'): -1.0,
              ref_model.metabolites.get_by_id('nadph_c'): 1.0,
              ref_model.metabolites.get_by_id('h_c'): 1.0,
              ref_model.metabolites.get_by_id('amet_c'): -1.0,
              ref_model.metabolites.get_by_id('ahcys_c'): 1.0,
              ref_model.metabolites.get_by_id('fad_c'): -1.0,
              ref_model.metabolites.get_by_id('fadh2_c'): 1.0,
              cofactor_metabolites_dict['mxcoa']: 1.0},

    'hpg_1': {ref_model.metabolites.get_by_id('pphn_c'): -1,
              ref_model.metabolites.get_by_id('34hpp_c'): 1,
              ref_model.metabolites.get_by_id('co2_c'): 1,
              ref_model.metabolites.get_by_id('h2o_c'): 1},

    'hpg_2': {ref_model.metabolites.get_by_id('34hpp_c'): -1,
              ref_model.metabolites.get_by_id('o2_c'): -1,
              ref_model.metabolites.get_by_id('h2o_c'): 1,
              ref_model.metabolites.get_by_id('4hmda_c'): 1},

    'hpg_3': {ref_model.metabolites.get_by_id('4hmda_c'): -1,
              ref_model.metabolites.get_by_id('fmn_c'): -1,
              ref_model.metabolites.get_by_id('fmnh2_c'): 1,
              ref_model.metabolites.get_by_id('nadh_c'): -1,
              ref_model.metabolites.get_by_id('nad_c'): 1,
              cofactor_metabolites_dict['4hbf']: 1},

    'hpg_4': {cofactor_metabolites_dict['4hbf']: -1,
              ref_model.metabolites.get_by_id('tyr__L_c'): -1,
              ref_model.metabolites.get_by_id('34hpp_c'): 1,
              cofactor_metabolites_dict['hpg']: 1},

    'bht': {ref_model.metabolites.get_by_id('tyr__L_c'): -1,
            ref_model.metabolites.get_by_id('o2_c'): -1,
            ref_model.metabolites.get_by_id('nadph_c'): -1,
            ref_model.metabolites.get_by_id('h_c'): -1,
            ref_model.metabolites.get_by_id('nadp_c'): 1,
            cofactor_metabolites_dict['hpg']: 1},

    'dpg': {ref_model.metabolites.get_by_id('accoa_c'): -1,
            ref_model.metabolites.get_by_id('malcoa_c'): -3,
            ref_model.metabolites.get_by_id('coa_c'): 4,
            ref_model.metabolites.get_by_id('co2_c'): 3,
            ref_model.metabolites.get_by_id('h2o_c'): 1,
            ref_model.metabolites.get_by_id('tyr__L_c'): -1,
            ref_model.metabolites.get_by_id('34hpp_c'): 1,
            cofactor_metabolites_dict['dpg']: 1},

    'gnat': {ref_model.metabolites.get_by_id('malcoa_c'): -1,
             ref_model.metabolites.get_by_id('co2_c'): 1,
             ref_model.metabolites.get_by_id('coa_c'): 1},

    'fkbh': {ref_model.metabolites.get_by_id('13dpg_c'): -1.0,
             ref_model.metabolites.get_by_id('pi_c'): 2.0},

    'ahba': {ref_model.metabolites.get_by_id('e4p_c'): -1.0,
             ref_model.metabolites.get_by_id('pep_c'): -1.0,
             ref_model.metabolites.get_by_id('pi_c'): 1.0,
             ref_model.metabolites.get_by_id('h2o_c'): 2.0},

    'acetyl': {ref_model.metabolites.get_by_id('accoa_c'): -1,
               ref_model.metabolites.get_by_id('coa_c'): 1,
               ref_model.metabolites.get_by_id('co2_c'): 1},

    'shikimic_acid': {ref_model.metabolites.get_by_id('skm_c'): -1,
                      ref_model.metabolites.get_by_id('co2_c'): 1},

    'fatty_acid': {ref_model.metabolites.get_by_id('accoa_c'): -1,
                   ref_model.metabolites.get_by_id('malcoa_c'): -3,
                   ref_model.metabolites.get_by_id('co2_c'): 4,
                   ref_model.metabolites.get_by_id('coa_c'): 4},

    'NH2': {ref_model.metabolites.get_by_id('malcoa_c'): -1,
            ref_model.metabolites.get_by_id('co2_c'): 2,
            ref_model.metabolites.get_by_id('gly_c'): -1,
            ref_model.metabolites.get_by_id('coa_c'): 1
            },

    'pip': {ref_model.metabolites.get_by_id('pyr_c'): -1,
            ref_model.metabolites.get_by_id('lys__L_c'): -1,
            ref_model.metabolites.get_by_id('h2o_c'): 1,
            ref_model.metabolites.get_by_id('nadph_c'): -1,
            ref_model.metabolites.get_by_id('nadp_c'): 1,
            ref_model.metabolites.get_by_id('h_c'): -1,
            cofactor_metabolites_dict['pip']: 1}

}


tailoring_metabolites_dict = {
    'glucose_6_phosphate': ref_model.metabolites.get_by_id('g6p_c'),
    '13biphosphoglycerate': ref_model.metabolites.get_by_id('13dpg_c'),
    'succinyl_coa': ref_model.metabolites.get_by_id('succoa_c'),
    'glycine': ref_model.metabolites.get_by_id('gly_c'),
    'coa': ref_model.metabolites.get_by_id('coa_c'),
    'co2': ref_model.metabolites.get_by_id('co2_c'),
    'atp': ref_model.metabolites.get_by_id('atp_c'),
    'amp': ref_model.metabolites.get_by_id('amp_c'),
    'ppi': ref_model.metabolites.get_by_id('ppi_c'),
    'pi': ref_model.metabolites.get_by_id('pi_c'),
    'nadh': ref_model.metabolites.get_by_id('nadph_c'),
    'nad': ref_model.metabolites.get_by_id('nadp_c'),
    'h+': ref_model.metabolites.get_by_id('h_c'),
    'h2o': ref_model.metabolites.get_by_id('h2o_c')
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

def _get_aa_formula(key):
    try:
        key_to_one_letter_dict = {value: key for key, value in one_letter_aa.items()}
        single_letter = key_to_one_letter_dict[key]
        return aa_formula[single_letter]
    except KeyError:
        return "na"