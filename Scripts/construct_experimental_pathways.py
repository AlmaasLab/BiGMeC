'''
Each function creates a set of metabolic reactions that together are able to synthesize the final metabolite

these reactions are returned as a list of cobra reactions.


'''


import cobra


# create metabolites that do not exist in ref_model
# but that are produced as intermediates in the reaction
cofactor_metabolites_dict = {
    'mxcoa': cobra.Metabolite('mxcoa_c', formula='C14H20N6O5S', name='Methoxymalonyl-CoA', compartment='c'),
    'bafA': cobra.Metabolite('bafilomycinA1_c', formula='C14H20N6O5S', name='Isobutyryl-CoA', compartment='c'),
    'fumamp': cobra.Metabolite('fumamp_c', formula='X', name='Fumaryl-AMP', compartment='c'),
    'mx__specific__acp_c':  cobra.Metabolite('mx__specific__acp_c', formula='X', name='methoxymalonate specific acyl carrier protein', compartment='c'),
    'c5n': cobra.Metabolite('c5n_c', formula='X', name='2-amino-3-hydroxycyclopent-2-enone', compartment='c'),
    'bafB': cobra.Metabolite('bafilomycinB1_c', formula='X', name='Bafilomycin B1', compartment='c'),
    'fumamp_pk': cobra.Metabolite('fumamp__pk_c', formula='X', 
                                  name='Fumaryl-AMP bound to polyketide', compartment='c'),
    'final_product': cobra.Metabolite('final_product_c', formula='X', name='Final product', compartment='c'),
    '3mbm': cobra.Metabolite('3mbmcoa_c', formula='C8H14O4', name='(3-Methylbutyl)malonic acid', compartment='c'),
    '5m2h': cobra.Metabolite('5m2hcoa__E_c', formula='X', name='	5-methyl-trans-hex-2-enoyl-ACP', compartment='c'),
    '5m3o': cobra.Metabolite('5mhcoa_c', formula='X', name='	5-Methyl-3-oxohexanoyl-ACP', compartment='c'),
    'i2b2': cobra.Metabolite('i2b2_c', formula='C8H12O4', name='	(2E)-2-Isobutyl-2-butenedioic acid', compartment='c'),
    'leupyrrin_1': cobra.Metabolite('leupyrrin_1_c', formula='X', name='leupyrrin intermediate', compartment='c'),
}


def create_baf_pathway(ref_model):
    reaction_metabolites = {
        ref_model.metabolites.get_by_id('nadph_c'): -11,
        ref_model.metabolites.get_by_id('nadp_c'): 11,
        ref_model.metabolites.get_by_id('h2o_c'): 4,
        ref_model.metabolites.get_by_id('coa_c'): 12,
        ref_model.metabolites.get_by_id('h_c'): -11,
        ref_model.metabolites.get_by_id('co2_c'): 12,
        cofactor_metabolites_dict['mxcoa']: -2,
        ref_model.metabolites.get_by_id('mmcoa__R_c'): -7,
        ref_model.metabolites.get_by_id('malcoa_c'): -2,
        ref_model.metabolites.get_by_id('ibcoa_c'): -1,
        cofactor_metabolites_dict['bafA']: 1
    }
    
    # Proposed biosynthetic pathway for FkbH-dependent (2R)-methoxymalonyl-ACP formation:
    mx_rxn_mets = {
        ref_model.metabolites.get_by_id('coa_c'): -1,
        ref_model.metabolites.get_by_id('13dpg_c'): -1,
        ref_model.metabolites.get_by_id('pi_c'): 2,
        ref_model.metabolites.get_by_id('nadp_c'): -1,
        ref_model.metabolites.get_by_id('nadph_c'): 1,
        ref_model.metabolites.get_by_id('h_c'): 1,
        ref_model.metabolites.get_by_id('amet_c'): -1,
        ref_model.metabolites.get_by_id('ahcys_c'): 1,
        ref_model.metabolites.get_by_id('fad_c'): -1,
        ref_model.metabolites.get_by_id('fadh2_c'): 1,
        cofactor_metabolites_dict['mxcoa']: 1
    }

    orf_3_fumamp_rxn_mets = {
        ref_model.metabolites.get_by_id('fum_c'): -1,
        ref_model.metabolites.get_by_id('atp_c'): -1,
        ref_model.metabolites.get_by_id('ppi_c'): 1,
        cofactor_metabolites_dict['fumamp']: 1,
    }

    baf_z_suc_gly_rxn_mets = {
        ref_model.metabolites.get_by_id('succoa_c'): -1,
        ref_model.metabolites.get_by_id('gly_c'): -1,
        ref_model.metabolites.get_by_id('co2_c'): 1,
        ref_model.metabolites.get_by_id('coa_c'): 1,
        ref_model.metabolites.get_by_id('5aop_c'): 1,
        ref_model.metabolites.get_by_id('h_c'): -1
    }

    orf_2_fumamp_pk_rxn_mets = {
        cofactor_metabolites_dict['bafA']: -1,
        cofactor_metabolites_dict['fumamp']: -1,
        ref_model.metabolites.get_by_id('amp_c'): 1,
        cofactor_metabolites_dict['fumamp_pk']: 1
    }

    baf_y_5aop_fumamp_pk_rxn_mets = {
        ref_model.metabolites.get_by_id('atp_c'): - 1,
        cofactor_metabolites_dict['fumamp_pk']: - 1,
        ref_model.metabolites.get_by_id('adp_c'): 1,
        ref_model.metabolites.get_by_id('pi_c'): 1,
        cofactor_metabolites_dict['bafB']: 1
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

    mx_rx = cobra.Reaction('mxcoa_synthesis')
    mx_rx.name = 'synthesis of methoxymalonyl-coa'
    mx_rx.lower_bound = 0.  # This is the default
    mx_rx.upper_bound = 1000.
    mx_rx.add_metabolites(mx_rxn_mets)

    ex_rx = cobra.Reaction('EX_final_product')
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.
    ex_rx.add_metabolites({cofactor_metabolites_dict['bafB']: -1})

    return [reaction, ex_rx, mx_rx, orf2_rx, orf3_rx, bafy_rx, bafz_rx]


def create_difficidin_pathway(ref_model):
    pk_metabolites = {
        ref_model.metabolites.get_by_id('nadph_c'): -14,
        ref_model.metabolites.get_by_id('nadp_c'): 14,
        ref_model.metabolites.get_by_id('h2o_c'): 7,
        ref_model.metabolites.get_by_id('coa_c'): 12,
        ref_model.metabolites.get_by_id('h_c'): -14,
        ref_model.metabolites.get_by_id('co2_c'): 12,
        ref_model.metabolites.get_by_id('malcoa_c'): -12,
        ref_model.metabolites.get_by_id('amet_c'): -3,
        ref_model.metabolites.get_by_id('ahcys_c'): 3,
        ref_model.metabolites.get_by_id('13dpg_c'): -1,
        ref_model.metabolites.get_by_id('pi_c'): 2,
        cofactor_metabolites_dict['final_product']: 1
    }

    pk_reaction = cobra.Reaction('Difficidin_synthesis')
    pk_reaction.name = 'Difficidin synthesis'
    pk_reaction.lower_bound = 0.  # This is the default
    pk_reaction.upper_bound = 1000.
    pk_reaction.add_metabolites(pk_metabolites)

    ex_rx = cobra.Reaction('EX_final_product')
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.
    ex_rx.add_metabolites({cofactor_metabolites_dict['final_product']: - 1})

    return [pk_reaction, ex_rx]


def create_anabaenopeptin_pathway(ref_model):
    pk_metabolites = {
        ref_model.metabolites.get_by_id('atp_c'): -6,
        ref_model.metabolites.get_by_id('amp_c'): 6,
        ref_model.metabolites.get_by_id('h2o_c'): 5,
        ref_model.metabolites.get_by_id('ppi_c'): 6,
        ref_model.metabolites.get_by_id('amet_c'): 1,
        ref_model.metabolites.get_by_id('ahcys_c'): 1,
        ref_model.metabolites.get_by_id('tyr__L_c'): -2,
        ref_model.metabolites.get_by_id('lys__L_c'): -1,
        ref_model.metabolites.get_by_id('ala__L_c'): -1,
        ref_model.metabolites.get_by_id('val__L_c'): -1,
        ref_model.metabolites.get_by_id('phe__L_c'): -1,
        cofactor_metabolites_dict['final_product']: 1
    }

    pk_reaction = cobra.Reaction('Anabaenopeptin_synthesis')
    pk_reaction.name = 'Anabaenopeptin synthesis'
    pk_reaction.lower_bound = 0.  # This is the default
    pk_reaction.upper_bound = 1000.
    pk_reaction.add_metabolites(pk_metabolites)

    ex_rx = cobra.Reaction('EX_final_product')
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.
    ex_rx.add_metabolites({cofactor_metabolites_dict['final_product']: - 1})

    return [pk_reaction, ex_rx]

def create_leupyrrin_pathway(ref_model):

    m3hhACP = {
        ref_model.metabolites.get_by_id('malcoa_c'): -1,
        ref_model.metabolites.get_by_id('ivcoa_c'): -1,
        ref_model.metabolites.get_by_id('co2_c'): 1,
        cofactor_metabolites_dict['5m3o']: 1
    }

    condensm3h = {
        cofactor_metabolites_dict['5m3o']: -1,
        ref_model.metabolites.get_by_id('h2o_c'): -1,
        cofactor_metabolites_dict['5m2h']: 1
    }
    # (3-Methylbutyl)malonic acid
    carbox = {
        ref_model.metabolites.get_by_id('nadph_c'): -1,
        ref_model.metabolites.get_by_id('h_c'): -1,
        ref_model.metabolites.get_by_id('nadp_c'): 1,
        ref_model.metabolites.get_by_id('co2_c'): -1,
        cofactor_metabolites_dict['3mbm']: 1,  # (3-Methylbutyl)malonic acid (basically final product.)
        cofactor_metabolites_dict['5m2h']: -1  #

    }
    '''These three are responsible for creating weird extender unit.'''

    pks_reaction = {
        ref_model.metabolites.get_by_id('pro__L_c'): -2,
        ref_model.metabolites.get_by_id('thr__L_c'): -1,
        ref_model.metabolites.get_by_id('malcoa_c'): -3,
        ref_model.metabolites.get_by_id('h2o_c'): 5,
        ref_model.metabolites.get_by_id('co2_c'): 4,
        ref_model.metabolites.get_by_id('coa_c'): 4,
        cofactor_metabolites_dict['3mbm']: -1,
        ref_model.metabolites.get_by_id('nadph_c'): -4,
        ref_model.metabolites.get_by_id('h_c'): -4,
        ref_model.metabolites.get_by_id('nadp_c'): 4,
        ref_model.metabolites.get_by_id('amet_c'): -2,
        ref_model.metabolites.get_by_id('ahcys_c'): 2,
        cofactor_metabolites_dict['leupyrrin_1']:1
    }

    otherrx_1 = {
        ref_model.metabolites.get_by_id('4mop_c'): -1,
        ref_model.metabolites.get_by_id('accoa_res_c'): -1,
        ref_model.metabolites.get_by_id('coa_c'): 1,
        cofactor_metabolites_dict['i2b2']: 1,
    }

    otherrx_2 = {
        cofactor_metabolites_dict['i2b2']: -1,
        cofactor_metabolites_dict['leupyrrin_1']: -1,
        cofactor_metabolites_dict['final_product']: 1
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

    otherrx_1_rx = cobra.Reaction('otherrx_1')
    otherrx_1_rx.lower_bound = 0.  # This is the default
    otherrx_1_rx.upper_bound = 1000.
    otherrx_1_rx.add_metabolites(otherrx_1)

    otherrx_2_rx = cobra.Reaction('otherrx_2')
    otherrx_2_rx.lower_bound = 0.  # This is the default
    otherrx_2_rx.upper_bound = 1000.
    otherrx_2_rx.add_metabolites(otherrx_2)

    ex_rx = cobra.Reaction('EX_final_product')
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.
    ex_rx.add_metabolites({cofactor_metabolites_dict['final_product']: - 1})

    return [m3hhACP_rx, condensm3h_rx, carbox_rx, pks_reaction_rx, otherrx_1_rx, otherrx_2_rx, ex_rx]


def create_tolaasin_pathway(ref_model):
    pk_metabolites = {
        ref_model.metabolites.get_by_id('thr__L_c'): -3,
        ref_model.metabolites.get_by_id('pro__L_c'): -1,
        ref_model.metabolites.get_by_id('ser__L_c'): -3,
        ref_model.metabolites.get_by_id('leu__L_c'): -4,
        ref_model.metabolites.get_by_id('val__L_c'): -4,
        ref_model.metabolites.get_by_id('gln__L_c'): -1,
        ref_model.metabolites.get_by_id('24dab_c'): -1,
        ref_model.metabolites.get_by_id('lys__L_c'): -1,
        ref_model.metabolites.get_by_id('atp_c'): -18,
        ref_model.metabolites.get_by_id('amp_c'): 18,
        ref_model.metabolites.get_by_id('ppi_c'): 18,
        cofactor_metabolites_dict['final_product']: 1
    }

    pk_reaction = cobra.Reaction('Tolaasin_synthesis')
    pk_reaction.name = 'Tolaasin synthesis'
    pk_reaction.lower_bound = 0.  # This is the default
    pk_reaction.upper_bound = 1000.
    pk_reaction.add_metabolites(pk_metabolites)

    ex_rx = cobra.Reaction('EX_final_product')
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.
    ex_rx.add_metabolites({cofactor_metabolites_dict['final_product']: - 1})

    return [pk_reaction, ex_rx]



def create_geldanamycin_pathway(ref_model):
    pk_metabolites = {
        ref_model.metabolites.get_by_id('nadph_c'): -14,
        ref_model.metabolites.get_by_id('nadp_c'): 14,
        ref_model.metabolites.get_by_id('h2o_c'): 7,
        ref_model.metabolites.get_by_id('coa_c'): 12,
        ref_model.metabolites.get_by_id('h_c'): -14,
        ref_model.metabolites.get_by_id('co2_c'): 12,
        ref_model.metabolites.get_by_id('malcoa_c'): -12,
        ref_model.metabolites.get_by_id('amet_c'): -3,
        ref_model.metabolites.get_by_id('ahcys_c'): 3,
        ref_model.metabolites.get_by_id('13dpg_c'): -1,
        ref_model.metabolites.get_by_id('pi_c'): 2,
        cofactor_metabolites_dict['final_product']: 1
    }

    pk_reaction = cobra.Reaction('Geldanamycin_synthesis')
    pk_reaction.name = 'Geldanamycin synthesis'
    pk_reaction.lower_bound = 0.  # This is the default
    pk_reaction.upper_bound = 1000.
    pk_reaction.add_metabolites(pk_metabolites)

    ex_rx = cobra.Reaction('EX_final_product')
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.
    ex_rx.add_metabolites({cofactor_metabolites_dict['final_product']: - 1})

    return [pk_reaction, ex_rx]

def create_oxazolo_pathway(ref_model):
    pk_metabolites = {
        ref_model.metabolites.get_by_id('nadph_c'): -8,
        ref_model.metabolites.get_by_id('ser__L_c'): -1,
        ref_model.metabolites.get_by_id('gly_c'): -2,
        ref_model.metabolites.get_by_id('nadp_c'): 8,
        ref_model.metabolites.get_by_id('h2o_c'): 8,
        ref_model.metabolites.get_by_id('coa_c'): 9,
        ref_model.metabolites.get_by_id('h_c'): -8,
        ref_model.metabolites.get_by_id('co2_c'): 10,
        ref_model.metabolites.get_by_id('malcoa_c'): -9,
        ref_model.metabolites.get_by_id('amet_c'): -5,
        ref_model.metabolites.get_by_id('ahcys_c'): 5,
        ref_model.metabolites.get_by_id('13dpg_c'): -1,
        ref_model.metabolites.get_by_id('pi_c'): 2,
        ref_model.metabolites.get_by_id('amet_c'): -1,
        ref_model.metabolites.get_by_id('ahcys_c'): 1,
        ref_model.metabolites.get_by_id('fad_c'): -1,
        ref_model.metabolites.get_by_id('fadh2_c'): 1,
        ref_model.metabolites.get_by_id('10fthf_c'): -1,
        ref_model.metabolites.get_by_id('thf_c'): -1,
        ref_model.metabolites.get_by_id('atp_c'): -3,
        ref_model.metabolites.get_by_id('amp_c'): 3,
        ref_model.metabolites.get_by_id('ppi_c'): 3,
        cofactor_metabolites_dict['final_product']: 1
    }

    pk_reaction = cobra.Reaction('Oxazolo_synthesis')
    pk_reaction.name = 'Oxazolo synthesis'
    pk_reaction.lower_bound = 0.  # This is the default
    pk_reaction.upper_bound = 1000.
    pk_reaction.add_metabolites(pk_metabolites)

    ex_rx = cobra.Reaction('EX_final_product')
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.
    ex_rx.add_metabolites({cofactor_metabolites_dict['final_product']: - 1})

    return [pk_reaction, ex_rx]

def create_oocydin_pathway(ref_model):

    pk_metabolites = {
        ref_model.metabolites.get_by_id('nadph_c'): -7, #ok
        ref_model.metabolites.get_by_id('nadp_c'): 7,#ok
        ref_model.metabolites.get_by_id('h2o_c'): 5,#ok
        ref_model.metabolites.get_by_id('coa_c'): 9,
        ref_model.metabolites.get_by_id('h_c'): -7,
        ref_model.metabolites.get_by_id('co2_c'): 9,
        ref_model.metabolites.get_by_id('malcoa_c'): -9,
        ref_model.metabolites.get_by_id('amet_c'): -2,
        ref_model.metabolites.get_by_id('ahcys_c'): 2,
        ref_model.metabolites.get_by_id('13dpg_c'): -1,
        ref_model.metabolites.get_by_id('pi_c'): 2,
        cofactor_metabolites_dict['final_product']: 1
    }

    pk_reaction = cobra.Reaction('Oocydin_synthesis')
    pk_reaction.name = 'Oocydin synthesis'
    pk_reaction.lower_bound = 0.  # This is the default
    pk_reaction.upper_bound = 1000.
    pk_reaction.add_metabolites(pk_metabolites)

    ex_rx = cobra.Reaction('EX_final_product')
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.
    ex_rx.add_metabolites({cofactor_metabolites_dict['final_product']: - 1})

    return [pk_reaction, ex_rx]

if __name__ == '__main__':
    ref_model_fn = "../Models/Sco-GEM.xml"
    model = cobra.io.read_sbml_model(ref_model_fn)
    create_oocydin_pathway(model)