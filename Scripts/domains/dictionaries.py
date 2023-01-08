#!/usr/bin/env python
# coding: utf-8
"""
Copyright 2020 Snorre Sulheim (snorre.sulheim@sintef.no)
https://github.com/AlmaasLab/BiGMeC

This file holds dictionaries used in the BiGMeC pipeline.

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
Updated: 08.01.2023
"""
import cobra
from pathlib import Path
# model_fn = '../Models/BiGG_universal_model.json'
# model_fn = '../Models/Sco-GEM.xml'



class BigmecDict(object):
    def __init__(self, model_fn):
        self._load_model(model_fn)
        self._make_dictionaries()

    def _load_model(self, model_fn):
        self.model_fn = model_fn
        if Path(model_fn).suffix == ".xml":
            self.model = cobra.io.read_sbml_model(model_fn)
        elif Path(model_fn).suffix == ".json":
            self.model = cobra.io.load_json_model(model_fn)
            if Path(model_fn).stem == "BiGG_universal_model":
                for m in self.model.metabolites:
                    m.compartment = "c"
        else:
            print("File type not recognized for reference model. Has to be json or sbml format.")
            raise IOError

    def get_met(self, m_id):
        return self.model.metabolites.get_by_id(m_id)

    def _raise_model_warning(model_fn):
        path = Path(model_fn)
        if model_fn.stem != "Sco-GEM.xml":
            warnings.warn("BiGMeC is not tested with tour selected reference model. \
                           This may cause unexpected errors. ")

    def _make_dictionaries(self):
        self.cofactor_metabolites_dict = {
            'nadph': self.model.metabolites.get_by_id('nadph_c'),
            'nadp': self.model.metabolites.get_by_id('nadp_c'),
            'h2o': self.model.metabolites.get_by_id('h2o_c'),
            'coa': self.model.metabolites.get_by_id('coa_c'),
            'co2': self.model.metabolites.get_by_id('co2_c'),
            'h+': self.model.metabolites.get_by_id('h_c'),
            'sam': self.model.metabolites.get_by_id('amet_c'),
            'sah': self.model.metabolites.get_by_id('ahcys_c'),
            'mxmal': cobra.Metabolite('mxmal_c', formula='C14H20N6O5S', name='Methoxymalonyl-ACP', compartment='c'),
            '4hbf': cobra.Metabolite('4hbf_c', formula='X', name='4-hydroxy-benzoyl-formate', compartment='c'),
            'hpg': cobra.Metabolite('4hpg_c', formula='X', name='4-hydroxy-phenyl-glycine', compartment='c'),
            'dpg': cobra.Metabolite('dpg_c', formula='X', name='dihydroxy-phenyl-glycine', compartment='c'),
            'bht': cobra.Metabolite('bht_c', formula='X', name='beta-hydroxy-tyrosine', compartment='c'),
            'pip': cobra.Metabolite('Lpipecol_c', formula='X', name='pipecolic acid', compartment='c'),
            'fatty_acid_X': cobra.Metabolite('fatty_acid_X_c', formula='X', 
                                        name='generic fatty acid for acylation in NRPS initiation', compartment='c'),
            'abu': cobra.Metabolite('2abu_c', formula='X', name='2-aminobutyrate', compartment='c'),
            'ahba': cobra.Metabolite('ahba_c', formula='C7H7NO3', name='3-Amino-5-hydroxybenzoate', compartment='c'), # https://www.genome.jp/dbget-bin/www_bget?C12107
            'bafA': cobra.Metabolite('bafilomycinA1_c', formula='X', name='Bafilomycin A1', compartment='c'),
            'bafB': cobra.Metabolite('bafilomycinB1_c', formula='X', name='Bafilomycin B1', compartment='c'),
            'fumamp': cobra.Metabolite('fumamp_c', formula='X', name='Fumaryl-AMP', compartment='c'),
            'mxacp':  cobra.Metabolite('mxacp_c', formula='X', name='methoxymalonate specific acyl carrier protein', compartment='c'),
            'c5n': cobra.Metabolite('c5n_c', formula='X', name='2-amino-3-hydroxycyclopent-2-enone', compartment='c'),
            'fumamp_pk': cobra.Metabolite('fumamp__pk_c', formula='X', 
                                          name='Fumaryl-AMP bound to polyketide', compartment='c'),
            'final_product': cobra.Metabolite('final_product_c', formula='X', name='Final product', compartment='c'),
            '3mbm': cobra.Metabolite('3mbmcoa_c', formula='C8H14O4', name='(3-Methylbutyl)malonic acid', 
                    compartment='c'),
            '5m2h': cobra.Metabolite('5m2hcoa__E_c', formula='X', name='5-methyl-trans-hex-2-enoyl-ACP', 
                    compartment='c'),
            '5m3o': cobra.Metabolite('5mhcoa_c', formula='X', name='5-Methyl-3-oxohexanoyl-ACP', compartment='c'),
            'i2b2': cobra.Metabolite('i2b2_c', formula='C8H12O4', name='(2E)-2-Isobutyl-2-butenedioic acid',
                                     compartment='c'),
            'leupyrrin_1': cobra.Metabolite('leupyrrin_1_c', formula='X', name='leupyrrin intermediate', 
                    compartment='c'),
            '4hmda': cobra.Metabolite('4hmda_c', name="(S)-4-Hydroxymandelate", compartment="c", charge="-1", formula="C8H7O4")}


        self.cofactor_reactions_dict = {  # reactions that are specific to certain domains
            'PKS_KR': {self.cofactor_metabolites_dict['nadph']: -1.0, self.cofactor_metabolites_dict['h+']: -1.0,
                       self.cofactor_metabolites_dict['nadp']: 1.0},
            'cMT': {self.cofactor_metabolites_dict['sam']: -1.0, self.cofactor_metabolites_dict['sah']: 1.0},
            'oMT': {self.cofactor_metabolites_dict['sam']: -1.0, self.cofactor_metabolites_dict['sah']: 1.0},
            'PKS_DH': {self.cofactor_metabolites_dict['h2o']: 1.0},
            'PKS_ER': {self.cofactor_metabolites_dict['nadph']: -1.0, self.cofactor_metabolites_dict['h+']: -1.0,
                       self.cofactor_metabolites_dict['nadp']: 1.0},
            'PKS_TE': {self.cofactor_metabolites_dict['h2o']: -1.0},
            'TD': {self.cofactor_metabolites_dict['nadph']: -1.0, self.cofactor_metabolites_dict['h+']: -1.0,
                       self.cofactor_metabolites_dict['nadp']: 1.0},
            'PKS_AT': {self.cofactor_metabolites_dict['coa']: 1.0, self.cofactor_metabolites_dict['co2']: 1.0},
            'Condensation': {self.cofactor_metabolites_dict['h2o']: -1.0},
            'nMT': {self.cofactor_metabolites_dict['sam']: -1.0, self.cofactor_metabolites_dict['sah']: 1.0},
            'mxmal': {
                      self.cofactor_metabolites_dict['mxacp']: -1.0,
                      self.model.metabolites.get_by_id('13dpg_c'): -1.0,
                      self.model.metabolites.get_by_id('pi_c'): 2.0,
                      self.model.metabolites.get_by_id('nad_c'): -1.0,
                      self.model.metabolites.get_by_id('nadh_c'): 1.0,
                      self.model.metabolites.get_by_id('h_c'): 1.0,
                      self.model.metabolites.get_by_id('amet_c'): -1.0,
                      self.model.metabolites.get_by_id('ahcys_c'): 1.0,
                      self.model.metabolites.get_by_id('fad_c'): -1.0,
                      self.model.metabolites.get_by_id('fadh2_c'): 1.0,
                      self.cofactor_metabolites_dict['mxmal']: 1.0},


            'hpg_1': {self.model.metabolites.get_by_id('pphn_c'): -1,
                      self.model.metabolites.get_by_id('34hpp_c'): 1,
                      self.model.metabolites.get_by_id('co2_c'): 1,
                      self.model.metabolites.get_by_id('h2o_c'): 1},

            'hpg_2': {self.model.metabolites.get_by_id('34hpp_c'): -1,
                      self.model.metabolites.get_by_id('o2_c'): -1,
                      self.model.metabolites.get_by_id('h2o_c'): 1,
                      self.cofactor_metabolites_dict['4hmda']: 1},

            'hpg_3': {self.cofactor_metabolites_dict['4hmda']: -1,
                      self.model.metabolites.get_by_id('fmn_c'): -1,
                      self.model.metabolites.get_by_id('fmnh2_c'): 1,
                      self.model.metabolites.get_by_id('nadh_c'): -1,
                      self.model.metabolites.get_by_id('nad_c'): 1,
                      self.cofactor_metabolites_dict['4hbf']: 1},

            'hpg_4': {self.cofactor_metabolites_dict['4hbf']: -1,
                      self.model.metabolites.get_by_id('tyr__L_c'): -1,
                      self.model.metabolites.get_by_id('34hpp_c'): 1,
                      self.cofactor_metabolites_dict['hpg']: 1},

            'bht': {self.model.metabolites.get_by_id('tyr__L_c'): -1,
                    self.model.metabolites.get_by_id('o2_c'): -1,
                    self.model.metabolites.get_by_id('nadph_c'): -1,
                    self.model.metabolites.get_by_id('h_c'): -1,
                    self.model.metabolites.get_by_id('nadp_c'): 1,
                    self.cofactor_metabolites_dict['hpg']: 1},

            'dpg': {self.model.metabolites.get_by_id('accoa_c'): -1,
                    self.model.metabolites.get_by_id('malcoa_c'): -3,
                    self.model.metabolites.get_by_id('coa_c'): 4,
                    self.model.metabolites.get_by_id('co2_c'): 3,
                    self.model.metabolites.get_by_id('h2o_c'): 1,
                    self.model.metabolites.get_by_id('tyr__L_c'): -1,
                    self.model.metabolites.get_by_id('34hpp_c'): 1,
                    self.cofactor_metabolites_dict['dpg']: 1},

            'gnat': {self.model.metabolites.get_by_id('malcoa_c'): -1,
                     self.model.metabolites.get_by_id('co2_c'): 1,
                     self.model.metabolites.get_by_id('coa_c'): 1},

            'fkbh': {self.model.metabolites.get_by_id('13dpg_c'): -1.0,
                     self.model.metabolites.get_by_id('pi_c'): 2.0},

            # See https://www.genome.jp/kegg-bin/show_pathway?rn01051 and
            # https://biocyc.org/META/NEW-IMAGE?type=PATHWAY&object=PWY-5979
            'ahba-synthesis': {self.model.metabolites.get_by_id('udpg_c'):-1,
                     self.model.metabolites.get_by_id('nad_c'): -1,
                     self.model.metabolites.get_by_id('nadh_c'): 1,
                     self.model.metabolites.get_by_id('h_c'): 3,
                     self.model.metabolites.get_by_id('gln__L_c'): -1,
                     self.model.metabolites.get_by_id('HC00591_c'): 1,
                     self.model.metabolites.get_by_id('udp_c'): 1,
                     self.model.metabolites.get_by_id('atp_c'): -1,
                     self.model.metabolites.get_by_id('adp_c'): 1,
                     self.model.metabolites.get_by_id('r5p_c'): -1,
                     self.model.metabolites.get_by_id('s7p_c') : 1,
                     self.model.metabolites.get_by_id('pep_c'): -1,
                     self.model.metabolites.get_by_id('pi_c'): 2,
                     self.cofactor_metabolites_dict["ahba"]: 1},

            'ahba': {self.cofactor_metabolites_dict["ahba"]: -1,
                     self.model.metabolites.get_by_id('atp_c'): -1,
                     self.model.metabolites.get_by_id('amp_c'):  1,
                     self.model.metabolites.get_by_id('ppi_c'):  1,
                     self.model.metabolites.get_by_id('h2o_c'):  1},

            'acetyl': {self.model.metabolites.get_by_id('accoa_c'): -1,
                       self.model.metabolites.get_by_id('coa_c'): 1,
                       self.model.metabolites.get_by_id('co2_c'): 1},

            'shikimic_acid': {self.model.metabolites.get_by_id('skm_c'): -1,
                              self.model.metabolites.get_by_id('atp_c'): -1,
                              self.model.metabolites.get_by_id('amp_c'):  1,
                              self.model.metabolites.get_by_id('ppi_c'):  1,
                              self.model.metabolites.get_by_id('h2o_c'):  1},

            'fatty_acid': {self.model.metabolites.get_by_id('accoa_c'): -1,
                           self.model.metabolites.get_by_id('malcoa_c'): -3,
                           self.model.metabolites.get_by_id('co2_c'): 4,
                           self.model.metabolites.get_by_id('coa_c'): 4},

            'NH2': {self.model.metabolites.get_by_id('malcoa_c'): -1,
                    self.model.metabolites.get_by_id('co2_c'): 2,
                    self.model.metabolites.get_by_id('gly_c'): -1,
                    self.model.metabolites.get_by_id('coa_c'): 1
                    },

            # https://www.sciencedirect.com/topics/immunology-and-microbiology/ascomycin
            'pip': {self.model.metabolites.get_by_id('pyr_c'): -1,
                    self.model.metabolites.get_by_id('lys__L_c'): -1,
                    self.model.metabolites.get_by_id('h2o_c'): 1,
                    self.model.metabolites.get_by_id('nadph_c'): -1,
                    self.model.metabolites.get_by_id('nadp_c'): 1,
                    self.model.metabolites.get_by_id('h_c'): -1,
                    self.cofactor_metabolites_dict['pip']: 1},

            'NRPS_acylating_loader': {self.cofactor_metabolites_dict["fatty_acid_X"]: -1,
                                      self.model.metabolites.get_by_id("coa_c"): 1},
            # R10992 in KEGG, could also be R10991
            'abu': {self.model.metabolites.get_by_id('2obut_c'): -1,
                    self.model.metabolites.get_by_id('ala__L_c'): -1,
                    self.model.metabolites.get_by_id('pyr_c'): 1,
                    self.cofactor_metabolites_dict["abu"]: 1},
        }

        # beta_hydroxy_acids = {
        #     '3hpp':  self.model.metabolites.get_by_id("3hpp_c"),
        #     'bhb': self.model.metabolites.get_by_id("bhb_c"),
        # }
        self.fatty_acyl_CoAs = {
            "3hpcoa": self.model.metabolites.get_by_id("3hpcoa_c"),
            # "ptpcoa": self.model.metabolites.get_by_id("ptpcoa_c"),
            "hxcoa":  self.model.metabolites.get_by_id("hxcoa_c"),
            # "hepcoa":  self.model.metabolites.get_by_id("hepcoa_c"),
            "occoa":  self.model.metabolites.get_by_id("occoa_c"),
            "dcacoa": self.model.metabolites.get_by_id("dcacoa_c")
        }

        # self.condensation_initiation_cofactors = {
        #     self.model.metabolites.get_by_id('h2o_c'): 1 
        # }

        self.tailoring_metabolites_dict = {
            'glucose_6_phosphate': self.model.metabolites.get_by_id('g6p_c'),
            '13biphosphoglycerate': self.model.metabolites.get_by_id('13dpg_c'),
            'succinyl_coa': self.model.metabolites.get_by_id('succoa_c'),
            'glycine': self.model.metabolites.get_by_id('gly_c'),
            'coa': self.model.metabolites.get_by_id('coa_c'),
            'co2': self.model.metabolites.get_by_id('co2_c'),
            'atp': self.model.metabolites.get_by_id('atp_c'),
            'amp': self.model.metabolites.get_by_id('amp_c'),
            'ppi': self.model.metabolites.get_by_id('ppi_c'),
            'pi': self.model.metabolites.get_by_id('pi_c'),
            'nadh': self.model.metabolites.get_by_id('nadph_c'),
            'nad': self.model.metabolites.get_by_id('nadp_c'),
            'h+': self.model.metabolites.get_by_id('h_c'),
            'h2o': self.model.metabolites.get_by_id('h2o_c')
        }

        self.tailoring_reactions_dict = {  # remember to add the secondary metabolite to these reactions
            'glycosyltransferase': {self.tailoring_metabolites_dict['glucose_6_phosphate']: -1,
                                    self.tailoring_metabolites_dict['h+']: 1,
                                    self.tailoring_metabolites_dict['pi']: 1
                                    },
            'glycerol': {self.tailoring_metabolites_dict['13biphosphoglycerate']: -1.0,
                         self.tailoring_metabolites_dict['pi']: 2.0,
                         self.tailoring_metabolites_dict['h+']: -1.0,
                         self.tailoring_metabolites_dict['nad']: 1.0,
                         self.tailoring_metabolites_dict['nadh']: -1.0,
                         },
            'ALA': {self.tailoring_metabolites_dict['succinyl_coa']: -1.0,
                    self.tailoring_metabolites_dict['glycine']: -1.0,
                    self.tailoring_metabolites_dict['atp']: -1.0,
                    self.tailoring_metabolites_dict['ppi']: 1.0,
                    self.tailoring_metabolites_dict['co2']: 1.0,
                    self.tailoring_metabolites_dict['coa']: 1.0,
                    self.tailoring_metabolites_dict['amp']: 1.0,
                    self.tailoring_metabolites_dict['h2o']: 1.0
                    }
                }

long_to_short = {'Malonyl-CoA': 'mal', 'Methylmalonyl-CoA': 'mmal', 'Methoxymalonyl-ACP': 'mxmal',
                 'Ethylmalonyl-CoA': 'emal', 'Isobutyryl-CoA': 'isobut', '2-Methylbutyryl-CoA': '2metbut',
                 'trans-1,2-CPDA': 'trans-1,2-CPDA', 'Acetyl-CoA': 'Acetyl-CoA', 'Benzoyl-CoA': 'benz',
                 'Propionyl-CoA': 'prop', '3-Methylbutyryl-CoA': '3metbut',
                 'CE-Malonyl-CoA': 'cemal', '2-Rhyd-Malonyl-CoA': '2Rhydmal', 'CHC-CoA': 'CHC-CoA',
                 'inactive': 'inactive'}

long_to_bigg = {'Malonyl-CoA': 'malcoa_c',
                'Methylmalonyl-CoA': 'mmcoa__R_c',
                'Methoxymalonyl-ACP': 'mxmal_c',
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
                'aoda': '99aoda_c',
                'UDP-glucose': 'udpg_c',
                'NAD+': 'nad_c',
                'NADH':  'nadh_c',
                'H+': 'h_c',
                'L-Glutamine': 'gln__L_c',
                '2-Oxoglutaramate': 'HC00591_c',
                'UDP': 'udp_c',
                'ATP': 'atp_c',
                'ADP': 'adp_c',
                'D-Ribose 5-phosphate': 'r5p_c',
                'Sedoheptulose 7-phosphate': 's7p_c',
                'Phosphoenolpyruvate': 'pep_c',
                'Orthophosphate': 'pi_c',
                 }

t1pks_extenders = {'Malonyl-CoA': 'C00083',
                   'Methylmalonyl-CoA': 'C00683',  # (S)-methylmalonyl-coa. (R) is C01213
                   'Methoxymalonyl-ACP': 'NOKEGG',
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