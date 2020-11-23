#!/usr/bin/env python
# coding: utf-8
"""
Copyright 2020 Snorre Sulheim (snorre.sulheim@sintef.no)
https://github.com/AlmaasLab/BiGMeC

This file is used to generate the manually curated pathways of the 8 BGCs used to evaluate the 
performace of BiGMeC.

This file is part of BiGMeC. BiGMeC is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. BiGMeC is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with BiGMeC.
If not, see <http://www.gnu.org/licenses/>.

Authors:
  - Snorre Sulheim, snorres.sulheim@sintef.no
  - Fredrik A. Fossheim

Date: 17.09.2020
"""

import cobra
from pathlib import Path
from dictionaries import *


def create_baf_pathway(ref_model):
    reaction_metabolites = {
        ref_model.metabolites.get_by_id('nadph_c'): -11,
        ref_model.metabolites.get_by_id('nadp_c'): 11,
        ref_model.metabolites.get_by_id('h2o_c'): 5,
        ref_model.metabolites.get_by_id('coa_c'): 10,
        ref_model.metabolites.get_by_id('h_c'): -11,
        ref_model.metabolites.get_by_id('co2_c'): 11,
        BDGLOBAL.cofactor_metabolites_dict['mxmal']: -2,
        BDGLOBAL.cofactor_metabolites_dict['mxacp']: 2,
        ref_model.metabolites.get_by_id('mmcoa__R_c'): -7,
        ref_model.metabolites.get_by_id('malcoa_c'): -2,
        ref_model.metabolites.get_by_id('ibcoa_c'): -1,
        BDGLOBAL.cofactor_metabolites_dict['bafA']: 1
    }
    

    orf_3_fumamp_rxn_mets = {
        ref_model.metabolites.get_by_id('fum_c'): -1,
        ref_model.metabolites.get_by_id('atp_c'): -1,
        ref_model.metabolites.get_by_id('ppi_c'): 1,
        BDGLOBAL.cofactor_metabolites_dict['fumamp']: 1,
    }

    baf_z_suc_gly_rxn_mets = {
        ref_model.metabolites.get_by_id('succoa_c'): -1,
        ref_model.metabolites.get_by_id('gly_c'): -1,
        ref_model.metabolites.get_by_id('co2_c'): 1,
        ref_model.metabolites.get_by_id('coa_c'): 1,
        ref_model.metabolites.get_by_id('5aop_c'): 1,
    }

    orf_2_fumamp_pk_rxn_mets = {
        BDGLOBAL.cofactor_metabolites_dict['bafA']: -1,
        BDGLOBAL.cofactor_metabolites_dict['fumamp']: -1,
        ref_model.metabolites.get_by_id('amp_c'): 1,
        BDGLOBAL.cofactor_metabolites_dict['fumamp_pk']: 1
    }

    baf_y_5aop_fumamp_pk_rxn_mets = {
        ref_model.metabolites.get_by_id('atp_c'): -1,
        BDGLOBAL.cofactor_metabolites_dict['fumamp_pk']: -1,
        ref_model.metabolites.get_by_id('5aop_c'): -1,
        ref_model.metabolites.get_by_id('adp_c'): 1,
        ref_model.metabolites.get_by_id('pi_c'): 1,
        BDGLOBAL.cofactor_metabolites_dict['bafB']: 1
    }

    orf3_rx = cobra.Reaction('ORF3')
    orf3_rx.lower_bound = 0.  # This is the default
    orf3_rx.upper_bound = 1000.
    orf3_rx.add_metabolites(orf_3_fumamp_rxn_mets)

    bafz_rx = cobra.Reaction('BAFZ')
    bafz_rx.lower_bound = 0.  # This is the default
    bafz_rx.upper_bound = 1000.
    bafz_rx.add_metabolites(baf_z_suc_gly_rxn_mets)

    orf2_rx = cobra.Reaction('ORF2')
    orf2_rx.lower_bound = 0.  # This is the default
    orf2_rx.upper_bound = 1000.
    orf2_rx.add_metabolites(orf_2_fumamp_pk_rxn_mets)

    bafy_rx = cobra.Reaction('BAFY')
    bafy_rx.lower_bound = 0.  # This is the default
    bafy_rx.upper_bound = 1000.
    bafy_rx.add_metabolites(baf_y_5aop_fumamp_pk_rxn_mets)

    reaction = cobra.Reaction('BAFILOMYCIN_SYNTHESIS')
    reaction.name = 'Bafilomycin synthesis'
    reaction.lower_bound = 0.  # This is the default
    reaction.upper_bound = 1000.
    reaction.add_metabolites(reaction_metabolites)

    mx_rx = cobra.Reaction('mxmal_synthesis')
    mx_rx.name = 'synthesis of methoxymalonyl-ACP'
    mx_rx.lower_bound = 0.  # This is the default
    mx_rx.upper_bound = 1000.
    mx_rx.add_metabolites(BDGLOBAL.cofactor_reactions_dict["mxmal"])

    ex_rx = cobra.Reaction('DM_secondary_metabolite')
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.
    ex_rx.add_metabolites({BDGLOBAL.cofactor_metabolites_dict['bafB']: -1})
    print(sum([reaction, ex_rx, orf2_rx, orf3_rx, bafy_rx, bafz_rx]))

    return [reaction, ex_rx, mx_rx, orf2_rx, orf3_rx, bafy_rx, bafz_rx]


def create_difficidin_pathway(ref_model):
    print("Difficidin")
    pk_metabolites = {
        ref_model.metabolites.get_by_id('nadph_c'): -14,
        ref_model.metabolites.get_by_id('nadp_c'): 14,
        ref_model.metabolites.get_by_id('h2o_c'): 5,
        ref_model.metabolites.get_by_id('coa_c'): 13,
        ref_model.metabolites.get_by_id('h_c'): -14,
        ref_model.metabolites.get_by_id('co2_c'): 12,
        ref_model.metabolites.get_by_id('malcoa_c'): -12,
        ref_model.metabolites.get_by_id('amet_c'): -3,
        ref_model.metabolites.get_by_id('ahcys_c'): 3,
        ref_model.metabolites.get_by_id('prpncoa_c'): -1,
        BDGLOBAL.cofactor_metabolites_dict['final_product']: 1
    }

    # prpncoa reaction
    #https://www.genome.jp/dbget-bin/www_bget?R12356
    prpncoa_reaction = cobra.Reaction("R12356")
    mets = {ref_model.metabolites.get_by_id('o2_c'): -1,
            ref_model.metabolites.get_by_id('h2o2_c'): 1,
            ref_model.metabolites.get_by_id('ppcoa_c'): -1,
            ref_model.metabolites.get_by_id('prpncoa_c'): 1,
            }
    prpncoa_reaction.add_metabolites(mets)
    prpncoa_reaction.bounds = (0, 1000)


    pk_reaction = cobra.Reaction('Difficidin_synthesis')
    pk_reaction.name = 'Difficidin synthesis'
    pk_reaction.lower_bound = 0.  # This is the default
    pk_reaction.upper_bound = 1000.
    pk_reaction.add_metabolites(pk_metabolites)

    ex_rx = cobra.Reaction('DM_secondary_metabolite')
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.
    ex_rx.add_metabolites({BDGLOBAL.cofactor_metabolites_dict['final_product']: - 1})
    print(sum([pk_reaction, ex_rx]))
    return [pk_reaction, ex_rx, prpncoa_reaction]


def create_anabaenopeptin_pathway(ref_model):
    pk_metabolites = {
        ref_model.metabolites.get_by_id('atp_c'): -6,
        ref_model.metabolites.get_by_id('amp_c'): 6,
        ref_model.metabolites.get_by_id('h2o_c'): 5,
        ref_model.metabolites.get_by_id('ppi_c'): 6,
        ref_model.metabolites.get_by_id('amet_c'): -1,
        ref_model.metabolites.get_by_id('ahcys_c'): 1,
        ref_model.metabolites.get_by_id('tyr__L_c'): -2,
        ref_model.metabolites.get_by_id('lys__L_c'): -1,
        ref_model.metabolites.get_by_id('ala__L_c'): -1,
        ref_model.metabolites.get_by_id('val__L_c'): -1,
        ref_model.metabolites.get_by_id('phe__L_c'): -1,
        BDGLOBAL.cofactor_metabolites_dict['final_product']: 1
    }

    pk_reaction = cobra.Reaction('Anabaenopeptin_synthesis')
    pk_reaction.name = 'Anabaenopeptin synthesis'
    pk_reaction.lower_bound = 0.  # This is the default
    pk_reaction.upper_bound = 1000.
    pk_reaction.add_metabolites(pk_metabolites)

    ex_rx = cobra.Reaction('DM_secondary_metabolite')
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.
    ex_rx.add_metabolites({BDGLOBAL.cofactor_metabolites_dict['final_product']: - 1})

    return [pk_reaction, ex_rx]

def create_leupyrrin_pathway(ref_model):

    m3hhACP = {
        ref_model.metabolites.get_by_id('malcoa_c'): -1,
        ref_model.metabolites.get_by_id('ivcoa_c'): -1,
        ref_model.metabolites.get_by_id('co2_c'): 1,
        ref_model.metabolites.get_by_id('coa_c'): 1,
        BDGLOBAL.cofactor_metabolites_dict['5m3o']: 1
    }

    condensm3h = {
        BDGLOBAL.cofactor_metabolites_dict['5m3o']: -1,
        ref_model.metabolites.get_by_id('h2o_c'): -1,
        BDGLOBAL.cofactor_metabolites_dict['5m2h']: 1
    }
    # (3-Methylbutyl)malonic acid
    carbox = {
        ref_model.metabolites.get_by_id('nadph_c'): -1,
        ref_model.metabolites.get_by_id('h_c'): -1,
        ref_model.metabolites.get_by_id('nadp_c'): 1,
        ref_model.metabolites.get_by_id('co2_c'): -1,
        BDGLOBAL.cofactor_metabolites_dict['3mbm']: 1,  # (3-Methylbutyl)malonic acid (basically final product.)
        BDGLOBAL.cofactor_metabolites_dict['5m2h']: -1  #

    }
    '''These three are responsible for creating weird extender unit.'''

    pks_reaction = {
        ref_model.metabolites.get_by_id('pro__L_c'): -2,
        ref_model.metabolites.get_by_id('fadh2_c'):  1,
        ref_model.metabolites.get_by_id('fad_c'): -1,
        ref_model.metabolites.get_by_id('pro__L_c'): -2,
        ref_model.metabolites.get_by_id('ppi_c'): 3,
        ref_model.metabolites.get_by_id('atp_c'): -3,
        ref_model.metabolites.get_by_id('amp_c'): 3,
        ref_model.metabolites.get_by_id('thr__L_c'): -1,
        ref_model.metabolites.get_by_id('malcoa_c'): -3,
        ref_model.metabolites.get_by_id('h2o_c'): 7,
        ref_model.metabolites.get_by_id('co2_c'): 4,
        ref_model.metabolites.get_by_id('coa_c'): 4,
        BDGLOBAL.cofactor_metabolites_dict['3mbm']: -1,
        ref_model.metabolites.get_by_id('nadph_c'): -4,
        ref_model.metabolites.get_by_id('h_c'): -4,
        ref_model.metabolites.get_by_id('nadp_c'): 4,
        ref_model.metabolites.get_by_id('amet_c'): -2,
        ref_model.metabolites.get_by_id('ahcys_c'): 2,
        BDGLOBAL.cofactor_metabolites_dict['leupyrrin_1']:1
    }

    # Skip leupyrrin tailoring
    # otherrx_1 = {
    #     ref_model.metabolites.get_by_id('4mop_c'): -1,
    #     ref_model.metabolites.get_by_id('accoa_c'): -1,
    #     ref_model.metabolites.get_by_id('coa_c'): 1,
    #     BDGLOBAL.cofactor_metabolites_dict['i2b2']: 1,
    # }

    otherrx_2 = {
        # BDGLOBAL.cofactor_metabolites_dict['i2b2']: -1,
        BDGLOBAL.cofactor_metabolites_dict['leupyrrin_1']: -1,
        BDGLOBAL.cofactor_metabolites_dict['final_product']: 1
    }

    m3hhACP_rx = cobra.Reaction('m3hhACP')
    m3hhACP_rx.lower_bound = 0.  # This is the default
    m3hhACP_rx.upper_bound = 1000.
    m3hhACP_rx.add_metabolites(m3hhACP)

    condensm3h_rx = cobra.Reaction('condensm3h')
    condensm3h_rx.lower_bound = 0.  # This is the default
    condensm3h_rx.upper_bound = 1000.
    condensm3h_rx.add_metabolites(condensm3h)

    carbox_rx = cobra.Reaction('carbox')
    carbox_rx.lower_bound = 0.  # This is the default
    carbox_rx.upper_bound = 1000.
    carbox_rx.add_metabolites(carbox)

    pks_reaction_rx = cobra.Reaction('pks_reaction')
    pks_reaction_rx.lower_bound = 0.  # This is the default
    pks_reaction_rx.upper_bound = 1000.
    pks_reaction_rx.add_metabolites(pks_reaction)

    # otherrx_1_rx = cobra.Reaction('otherrx_1')
    # otherrx_1_rx.lower_bound = 0.  # This is the default
    # otherrx_1_rx.upper_bound = 1000.
    # otherrx_1_rx.add_metabolites(otherrx_1)

    otherrx_2_rx = cobra.Reaction('otherrx_2')
    otherrx_2_rx.lower_bound = 0.  # This is the default
    otherrx_2_rx.upper_bound = 1000.
    otherrx_2_rx.add_metabolites(otherrx_2)

    ex_rx = cobra.Reaction('DM_secondary_metabolite')
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.
    ex_rx.add_metabolites({BDGLOBAL.cofactor_metabolites_dict['final_product']: - 1})
    print("Leupyrrin")
    print(sum([pks_reaction_rx, otherrx_2_rx, ex_rx]))
    return [m3hhACP_rx, condensm3h_rx, carbox_rx, pks_reaction_rx,  otherrx_2_rx, ex_rx] #otherrx_1_rx,


def create_tolaasin_pathway(ref_model):
    # Necessary to enable consumption of 3hocoa_c
    hocoa_mets = {ref_model.metabolites.get_by_id('3oocoa_c'): -1.0,
                  ref_model.metabolites.get_by_id('coa_c'): -1.0,
                  ref_model.metabolites.get_by_id('accoa_c'): 1.0,
                  ref_model.metabolites.get_by_id('hxcoa_c'): 1.0}

    hocoa_rx = cobra.Reaction('hocoa_pathway')
    hocoa_rx.lower_bound = -1000. 
    hocoa_rx.upper_bound = 1000.
    hocoa_rx.add_metabolites(hocoa_mets)


    # Incoporarate threonine and not 2,3-didehydrobutyrine, see https://doi.org/10.1186/s12934-016-0502-y
    pk_metabolites = {
        ref_model.metabolites.get_by_id('coa_c'): 1,
        ref_model.metabolites.get_by_id('3hocoa_c'): -1,
        ref_model.metabolites.get_by_id('thr__L_c'): -3,
        ref_model.metabolites.get_by_id('pro__L_c'): -1,
        ref_model.metabolites.get_by_id('ser__L_c'): -2,
        ref_model.metabolites.get_by_id('hom__L_c'): -1,
        ref_model.metabolites.get_by_id('leu__L_c'): -4,
        ref_model.metabolites.get_by_id('val__L_c'): -4,
        ref_model.metabolites.get_by_id('gln__L_c'): -1,
        ref_model.metabolites.get_by_id('24dab_c'): -1,
        ref_model.metabolites.get_by_id('lys__L_c'): -1,
        ref_model.metabolites.get_by_id('atp_c'): -18,
        ref_model.metabolites.get_by_id('amp_c'): 18,
        ref_model.metabolites.get_by_id('ppi_c'): 18,
        ref_model.metabolites.get_by_id('h2o_c'): 17,
        BDGLOBAL.cofactor_metabolites_dict['final_product']: 1
    }

    pk_reaction = cobra.Reaction('Tolaasin_synthesis')
    pk_reaction.name = 'Tolaasin synthesis'
    pk_reaction.lower_bound = 0.  # This is the default
    pk_reaction.upper_bound = 1000.
    pk_reaction.add_metabolites(pk_metabolites)

    ex_rx = cobra.Reaction('DM_secondary_metabolite')
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.
    ex_rx.add_metabolites({BDGLOBAL.cofactor_metabolites_dict['final_product']: - 1})
    print(sum([pk_reaction, ex_rx]))
    return [pk_reaction, ex_rx, hocoa_rx]



def create_geldanamycin_pathway(ref_model):
    pk_metabolites = {
        BDGLOBAL.cofactor_metabolites_dict['mxmal']: -2,
        BDGLOBAL.cofactor_metabolites_dict['mxacp']: 2,
        BDGLOBAL.cofactor_metabolites_dict['ahba']: -1,
        ref_model.metabolites.get_by_id('atp_c'): -1,
        ref_model.metabolites.get_by_id('amp_c'):  1,
        ref_model.metabolites.get_by_id('mmcoa__R_c'): -4,
        ref_model.metabolites.get_by_id('nadph_c'): -10,
        ref_model.metabolites.get_by_id('nadp_c'): 10,
        ref_model.metabolites.get_by_id('h2o_c'): 6,
        ref_model.metabolites.get_by_id('coa_c'): 5,
        ref_model.metabolites.get_by_id('h_c'): -10,
        ref_model.metabolites.get_by_id('co2_c'): 7,
        ref_model.metabolites.get_by_id('malcoa_c'): -1,
        ref_model.metabolites.get_by_id('ppi_c'): 1,
        BDGLOBAL.cofactor_metabolites_dict['final_product']: 1
    }
    
    
    mx_rx = cobra.Reaction('mxmal_synthesis')
    mx_rx.name = 'synthesis of methoxymalonyl-ACP'
    mx_rx.lower_bound = 0.  # This is the default
    mx_rx.upper_bound = 1000.
    mx_rx.add_metabolites(BDGLOBAL.cofactor_reactions_dict["mxmal"])

    ahba_synthesis_reaction = cobra.Reaction('ahba_synthesis')
    ahba_synthesis_reaction.name = "AHBA synthesis"
    ahba_synthesis_reaction.bounds = (0,1000)
    ahba_synthesis_reaction.add_metabolites(BDGLOBAL.cofactor_reactions_dict['ahba-synthesis'])
    

    pk_reaction = cobra.Reaction('Geldanamycin_synthesis')
    pk_reaction.name = 'Geldanamycin synthesis'
    pk_reaction.lower_bound = 0.  # This is the default
    pk_reaction.upper_bound = 1000.
    pk_reaction.add_metabolites(pk_metabolites)

    ex_rx = cobra.Reaction('DM_secondary_metabolite')
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.
    ex_rx.add_metabolites({BDGLOBAL.cofactor_metabolites_dict['final_product']: - 1})
    print(sum([pk_reaction, ex_rx]))
    return [mx_rx, pk_reaction, ex_rx, ahba_synthesis_reaction]

def create_oxazolo_pathway(ref_model):
    pk_metabolites = {
        ref_model.metabolites.get_by_id('for_c'): -1,
        ref_model.metabolites.get_by_id('gly_c'): -2,
        ref_model.metabolites.get_by_id('h2o_c'): 11,
        ref_model.metabolites.get_by_id('atp_c'): -3,
        ref_model.metabolites.get_by_id('amp_c'): 3,
        ref_model.metabolites.get_by_id('ppi_c'): 3,
        ref_model.metabolites.get_by_id('malcoa_c'): -9,
        ref_model.metabolites.get_by_id('co2_c'): 10,
        ref_model.metabolites.get_by_id('coa_c'): 9,
        ref_model.metabolites.get_by_id('nadph_c'): -9,
        ref_model.metabolites.get_by_id('h_c'): -9,
        ref_model.metabolites.get_by_id('nadp_c'): 9,
        ref_model.metabolites.get_by_id('amet_c'): -5,
        ref_model.metabolites.get_by_id('ahcys_c'): 5,
        BDGLOBAL.cofactor_metabolites_dict['mxmal']: -1,
        BDGLOBAL.cofactor_metabolites_dict['mxacp']: 1,
        ref_model.metabolites.get_by_id('ser__L_c'): -1,
        BDGLOBAL.cofactor_metabolites_dict['final_product']: 1
    }
    mx_rx = cobra.Reaction('mxmal_synthesis')
    mx_rx.name = 'synthesis of methoxymalonyl-ACP'
    mx_rx.lower_bound = 0.  # This is the default
    mx_rx.upper_bound = 1000.
    mx_rx.add_metabolites(BDGLOBAL.cofactor_reactions_dict["mxmal"])

    pk_reaction = cobra.Reaction('Oxazolo_synthesis')
    pk_reaction.name = 'Oxazolo synthesis'
    pk_reaction.lower_bound = 0.  # This is the default
    pk_reaction.upper_bound = 1000.
    pk_reaction.add_metabolites(pk_metabolites)

    ex_rx = cobra.Reaction('DM_secondary_metabolite')
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.
    ex_rx.add_metabolites({BDGLOBAL.cofactor_metabolites_dict['final_product']: - 1})

    return [mx_rx, pk_reaction, ex_rx]

def create_oocydin_pathway(ref_model):

    pk_metabolites = {
        ref_model.metabolites.get_by_id('nadph_c'): -7, #ok
        ref_model.metabolites.get_by_id('nadp_c'): 7,#ok
        ref_model.metabolites.get_by_id('h2o_c'): 6,#ok
        ref_model.metabolites.get_by_id('coa_c'): 9,
        ref_model.metabolites.get_by_id('h_c'): -7,
        ref_model.metabolites.get_by_id('co2_c'): 9,
        ref_model.metabolites.get_by_id('malcoa_c'): -9,
        ref_model.metabolites.get_by_id('amet_c'): -2,
        ref_model.metabolites.get_by_id('ahcys_c'): 2,
        ref_model.metabolites.get_by_id('13dpg_c'): -1,
        ref_model.metabolites.get_by_id('pi_c'): 2,
        BDGLOBAL.cofactor_metabolites_dict['final_product']: 1
    }

    pk_reaction = cobra.Reaction('Oocydin_synthesis')
    pk_reaction.name = 'Oocydin synthesis'
    pk_reaction.lower_bound = 0.  # This is the default
    pk_reaction.upper_bound = 1000.
    pk_reaction.add_metabolites(pk_metabolites)

    ex_rx = cobra.Reaction('DM_secondary_metabolite')
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.
    ex_rx.add_metabolites({BDGLOBAL.cofactor_metabolites_dict['final_product']: - 1})
    print(sum([pk_reaction, ex_rx]))
    return [pk_reaction, ex_rx]

def make_json_pathways(ref_model):

    pathway_dict = {
        "Bafilomycin": create_baf_pathway,
        "Difficidin": create_difficidin_pathway,
        "Anabaenopeptin": create_anabaenopeptin_pathway,
        "Tolaasin": create_tolaasin_pathway,
        "Leupyrrin": create_leupyrrin_pathway,
        "Geldanamycin": create_geldanamycin_pathway,
        "Oocydin": create_oocydin_pathway,
        "Oxazolomycin": create_oxazolo_pathway
    }
    folder = Path("../Data/validation_pathways")
    folder.mkdir(exist_ok=True)
    
    for key, value in pathway_dict.items():
        fn = folder / "{0}.json".format(key)
        model = cobra.Model()
        reactions = value(ref_model)
        model.add_reactions(reactions)
        cobra.io.save_json_model(model, str(fn))



if __name__ == '__main__':
    ref_model_fn = "../Models/Sco-GEM.xml"

    global BDGLOBAL
    BDGLOBAL = BigmecDict(ref_model_fn)
    ref_model = BDGLOBAL.model
    if 0:
        create_oocydin_pathway(ref_model)
    if 1:
        make_json_pathways(ref_model)