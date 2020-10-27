#!/usr/bin/env python
# coding: utf-8
"""
Authors:
  - Snorre Sulheim, snorres.sulheim@sintef.no
  - Fredrik Fossheim

Date: 17.09.2020
Lisence: CC-BY-4.0

This is the main file used to run BiGMeC to predict metabolic pathways from identified and annotated BGCs.

"""

from Bio import SeqIO
import re
import os
import json
import sys
import re
import pandas as pd
import cobra  
import copy
from itertools import groupby
import warnings
from pathlib import Path
import numpy as np

from domains import *
from dictionaries import *


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
                return {'start': low, 'end': high,
                         'type': typer,
                         'core_number': core['core_number']}
            elif core['type'] in RiPPs:
                return core

    return core


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
    # Get the range of all core genes
    low = min([core["start"] for core in core_list])
    high = max([core["end"] for core in core_list])

    merged_core = merge_cores(core_list, high, low)

    return merged_core


def find_cores_in_cluster(gb_list):
    '''
    In some files, the CDS shows up before the proto core, so tou actually have to find proto_core on beforehand
    (this is in order to be able to know which genes are part of the PKS synthesis)
    This is the reason that we iterate through the gb_list one time before the big next iteration
    :param gb_list:
    :return: A list of cores (cores are overlapping. we are simply interested in the core genes for this project )
    '''
    core_list = []
    
    # backup_core_list is there because some BGCs do not have any locations for their clusters. 
    # this i think is just for older BGCs that are in mibig, as this does not always happen.
    backup_core_list = []
    
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


def structure_gbk_information(merged_core, gb_list):
    CDS = []
    core_domains = {}
    CDS_number = -1
    for gb_record in gb_list:
        for feat in gb_record.features:
            # print(feat, CDS_number, feat.qualifiers)
            if feat.type == 'CDS':
                if feat.qualifiers.get('gene_functions'):
                    smcog = []
                    gene_function = feat.qualifiers.get('gene_functions')
                    for single in gene_function:
                        if 'SMCOG' in single:
                            smcog.append(single.split('SMCOG')[1].split(':')[0])

                else:
                    smcog = ['No_smCOG']

                cds_is_core = False
                try:
                    if feat.location.start >= merged_core['start'] and feat.location.end <= merged_core['end']:
                        cds_is_core = True
                except AttributeError:
                    warnings.warn('Core has no location. Unsure of consequences')
                    # Instantiate CDS-object, consisting of:
                    # list of domains, gene_ontologies, smCOG, EC numbers (possible with many because of many GOs, rxn_in,
                    # rxn_out, core_gene):


                try:
                    CDS.append({'smcog': smcog, 'core_gene': cds_is_core, 
                                'domains': [], 'strand': feat.location.strand})
                    CDS_number += 1
                except AttributeError:
                    warnings.warn('gene had no location')
                
            elif feat.type == 'aSModule':
                if feat.qualifiers.get('monomer_pairings'):
                    aids = feat.qualifiers.get('monomer_pairings')
                    #for core in core_list:
                    if feat.location.start >= merged_core['start'] and feat.location.end <= merged_core['end']:
                        if merged_core['core_number'] not in core_domains:
                            core_domains[merged_core['core_number']] = {'type': merged_core['type'], 'modules': [
                                {'extender_unit': copy.copy(aids)[0].replace('&gt;', '>'),
                                 'start': feat.location.start,
                                 'end': feat.location.end}]}
                            #  replace('&gt;', '>') -> sometimes there is some parsing error that cant read '>'
                        else:
                            core_domains[merged_core['core_number']]['modules'] += [
                                {'extender_unit': copy.copy(aids)[0].replace('&gt;', '>'),
                                 'start': feat.location.start,
                                 'end': feat.location.end}]

            elif feat.type == 'CDS_motif':
                if feat.qualifiers.get('core_sequence'):
                    core_domains[len(core_domains) + 1] = {'type': feat.qualifiers.get('peptide')[0],
                                                           'RiPP': feat.qualifiers.get('core_sequence')[0]}
            elif feat.type == 'aSDomain':                                   
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

                elif feat.qualifiers.get('aSDomain') == ['PKS_KR']:  # KR domains
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
                                {'type': 'CAL_domain',
                                 'activity': CAL_prediction.split('Minowa: ')[1],
                                 'start': feat.location.start,
                                 'end': feat.location.end,
                                 'gene': feat.qualifiers.get('locus_tag'),
                                 'core_gene': cds_is_core}
                            )
                elif feat.qualifiers.get('aSDomain') == ["AMP-binding"]:
                    specificity = feat.qualifiers.get('specificity')[0].split(":")[-1].strip()
                    CDS[CDS_number]['domains'].append(
                        {'type': "AMP-binding",
                         'activity': True,
                         'start': feat.location.start,
                         'end': feat.location.end,
                         'gene': feat.qualifiers.get('locus_tag'),
                         'core_gene': cds_is_core,
                         'AT_specificity': specificity}
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

def _get_cluster_type(gb_list, merged_core):
    types = []
    for gb in gb_list:
        for f in gb.features:
            if f.type == "region":
                types += f.qualifiers["product"]

    return list(set(types))


def parse_antismash_gbk(gbk_path, json_path, save_json = True):
    gb_list = get_gb_list_from_antismash_output(gbk_path)
    core_list = find_cores_in_cluster(gb_list)
    merged_core = merge_core_list(core_list)
    cluster_type = _get_cluster_type(gb_list, merged_core)
    CDS = structure_gbk_information(merged_core, gb_list)
    CDS["Cluster types"] = cluster_type

    if save_json:
        with open(json_path, 'w+') as f:
            json.dump(CDS, f)
        f.close()
    return CDS


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


def fix_transat_modules(domain_or_module):
    '''
    Removes extra acp domains (?)

    :param domain_or_module: a list of domains and modules that antismash have found. This is the first of several
    "fixings" of this and only applies to transAT BGCs.
    
    :return: domain_or_module, except some modules have been broken down if they are ks-
    '''
    
    N = len(domain_or_module)
    for i in range(N - 1):
        module = domain_or_module[i]
        domain = domain_or_module[i + 1]
        if module['element'] == 'module' and domain['element'] == 'domain':
            if domain['info']['type'] in reducing_domains:
                N_d = len(module['domains'])
                for j in range(N_d - 1):
                    # insert domains at end, pop at start
                    if module['domains'][j]['type'] in acp_domains:
                        if module['domains'][j+1]['type'] in at_domains:

                            pre_domains  = [{'element': 'domain', 'info': x} for x in module['domains'][:j]]
                            post_domains = [{'element': 'domain', 'info': x} for x in module['domains'][j + 1:]]
                            domain_or_module = fix_transat_modules(domain_or_module[0:i] + pre_domains + 
                                                post_domains + domain_or_module[i + 1:])
                return domain_or_module
    return domain_or_module


def find_and_replace_DHD_domains(domain_or_module):
    '''
    :param domain_or_module: sequence of domains and modules
    :return: find all domain sequences that are KS-DH-ACP and span across two genes. 
    These we put into a non-extending module with extender_unit: 'DHD'
    '''
    DHD_modules = []
    for index in range(len(domain_or_module)-2):
        three_domains = [domain_or_module[index], domain_or_module[index+1], domain_or_module[index+2]]
        # Check if all three are domains
        if np.all([d["element"] == "domain" for d in three_domains]):
            d_types = [d["info"]["type"] for d in three_domains]
            if (d_types[0] == 'PKS_KS') and (d_types[1] in dh_domains) and (d_types[2] in acp_domains):
                # check that KS and DH is actually on two different genes
                if domain_or_module[index]['info']['gene'] != domain_or_module[index + 1]['info']['gene']:
                    DHD_module = {'element': 'module',
                                  'info': {'extender_unit': 'DHD', 
                                           'start': domain_or_module[index]['info']['start'],
                                           'end':   domain_or_module[index+2]['info']['end']},
                                  'domains': [d['info'] for d in three_domains]}
                    DHD_modules.append([index, DHD_module])              
                    # domain_or_module = domain_or_module[0:index] + [new_module] + domain_or_module[index + 3:]
    if len(DHD_modules):
        all_domains = []
        prev_idx = 0
        for index, new_module in DHD_modules:
            if index < prev_idx:
                raise IndexError
            all_domains += domain_or_module[prev_idx:index]
            all_domains.append(new_module)
            prev_idx = index + 3
        all_domains += domain_or_module[prev_idx:]
    else:
        all_domains = domain_or_module
    return all_domains


def deactivate_dh_domains_in_DHD_modules(domain_or_module):
    """
    The DH domains in DHD modules seems to not be active, and they are therefore deactivated
    """
    for module in domain_or_module:
        if module['element'] == 'domain':
            # Only interested in modules
            continue

        if module["info"]["extender_unit"] == "DHD":
            for domain in module["domains"]:
                if domain["type"] in dh_domains:
                    domain["activity"] = False
    return domain_or_module



def deactivate_reducing_domains_when_kr_is_inactive(domain_or_module, deactivate_if_no_KR):
    '''
    If there is a KR domain in a module the activity of other reducing domains is set to the same activity (True or False). This is because the activity is only predicted for the KR domain, but if the KR domain is inactive, the other reducing domains are neither active (They need to act on the OH-group created by the KR domain).

    :param domain_or_module
        - List of domains and modules
    :param deactivate_if_no_KR
        - bool, determine the action if no KR domain is present in the module. In T1PKS we would like to set this to True, but for trans-AT we set it false since it is possible that a KR domain acts in trans. 

    :return: List of domains and module with corrected activity
    '''
    for module in domain_or_module:
        if module['element'] == 'domain':
            # Only interested in modules
            continue

        # Find KR domain
        KR_domain = False
        KR_activity = False
        for domain in module["domains"]:
            if domain['type'] == 'PKS_KR':
                KR_domain = True
                KR_activity = domain['activity']

        # Skip to next module if there is no KR domain in this module and we don't want to deactive other domains 
        # in case of missing KR domain 
        if KR_domain or deactivate_if_no_KR:
            # Set activity of other DH / ER domains to the activity of the KR domain
             for domain in module["domains"]:
                if domain['type'] in dh_er_domains:
                    domain['activity'] = KR_activity
        else:
            continue


    return domain_or_module


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
                new_domains = [domain_or_module[i]['info'] for i in range(start, end)]
                new_module = {'element': 'module', 
                              'info': {'extender_unit': 'mal', 
                                       'start': res[start]['info']['start'],
                                       'end': res[end]['info']['end']},
                              'domains': new_domains}
                all_modules = res[0:start] + [new_module] + res[end + 1:]
                res = find_and_replace_cut_transat_modules(all_modules)
                return res
    return res


def add_rogue_transat_domains(module):
    """
    Removes duplicate acp-domains from a domain or module.
    This phenoenon is commonly observed for transAT-PKS and from anecdotal evidence they never have any effect on
    the biosynthesis of the secondary metabolite. The reason they are removed is that this makes it easier to
    look for modules that are wrongly predicted by antismash
    """

    for ind, domain in enumerate(module):
        if domain['element'] == 'domain':
            if domain['info']['type'] in acp_domains:
                module.pop(ind)

    for index in range(1, len(module) - 1):
        if module[index]['element'] == 'module' and module[index + 1]['element'] == 'domain':

            if module[index + 1]['info']['type'] in reducing_domains and not module[index]['info']['extender_unit'] == 'DHD':
                # case where domain is located right next to module
                not_domain_exists_in_prev_mod = True
                for donk in module[index]['domains']:
                    if module[index + 1]['info']['type'] == donk['type']:
                        # check to see if the rogue domain exists in the module we want to add it to
                        not_domain_exists_in_prev_mod = False
                if not_domain_exists_in_prev_mod:
                    module[index]['domains'].append(module[index + 1]['info'])

            elif module[index + 1]['info']['type'] in reducing_domains and module[index]['info']['extender_unit'] == 'DHD' and \
                    module[index - 1]['element'] == 'module':
                # case where domain is located right next to DHD module, located right next to another module
                not_domain_exists_in_prev_mod = True
                for donk in module[index - 1]['domains']:
                    if module[index + 1]['info']['type'] == donk['type']:
                        # check to see if the rogue domain exists in the module we want to add it to
                        not_domain_exists_in_prev_mod = False

                if not_domain_exists_in_prev_mod:
                    module[index - 1]['domains'].append(module[index + 1]['info'])

    for ind, a in enumerate(module):
        if ind < len(module) - 1:
            if a['element'] == 'module':
                for dom in a['domains']:
                    if dom == module[ind + 1]['info']:
                        module.pop(ind + 1)
                    elif ind < len(module) - 2:
                        if dom == module[ind + 2]['info']:
                            module.pop(ind + 2)

    return module


def add_transat_metabolic_pathway(data, model, core_number, tailoring_reactions):
    """
    :param data: all data that has been structured from gbk-file
    :return: modules max: The maximum detected modules, following rule that it needs a KS-domain and an ACP domain
    Modules mid: All modules that contain KS, AT and ACP domain.
    the modules_min is
    """

    
    # domains x modules is now a dict sorted by placement on BGC.
    # do this to find out if domains are found in this sequence: ks-at-acp
    domains, actual_domains, domains_x_modules = _get_domains(data, core_number)

    # assume that DHD-domains are inactive, however, may still dehydratase the previous domain(?)
    # assue that all others are active.
    domains_x_modules = fix_transat_modules(domains_x_modules)
    domains_x_modules = find_and_replace_DHD_domains(domains_x_modules)
    domains_x_modules = find_and_replace_load_modules(domains_x_modules)
    
    domains_x_modules = find_and_replace_cut_transat_modules(domains_x_modules)
    
    domains_x_modules = add_rogue_transat_domains(domains_x_modules)

    # KR domains are basically always active
    # domains_x_modules = activate_kr_domains(domains_x_modules)

    # Set the activity of the DH / ER domains 
    domains_x_modules = deactivate_reducing_domains_when_kr_is_inactive(domains_x_modules, False)
    domains_x_modules = deactivate_dh_domains_in_DHD_modules(domains_x_modules)
    # Remove omt modules. 
    # these "rogue" omt domains are usually for tailoring reactions of the pks and is therefore just a bonus
    # It is a weird case, but should make it more correct for the most part, although the methodology behind
    # it is completely wrong.
    domains_x_modules = search_and_destroy_omt_modules(domains_x_modules)
    model = create_t1_transat_nrps_model(data['core_structure'][core_number], domains_x_modules, model,
                                         tailoring_reactions)

    return model



def _is_standard_load_module(module):
    """
    Determines if the module is a classical load module with the A/AT-ACP/PCP-domains
    """
    
    if module['element'] == 'module':
        if len(module['domains']) == 2:
            if module['domains'][0]['type'] in loader_at_domains:
                if module['domains'][1]['type'] in loader_acp_domains:
                    return True
    return False

def find_and_replace_load_modules(domain_or_module):
    '''
    :param domain_or_module: list of modules and domains
    :return: domain or module, but add a starter module.
    
    Function that identifies / add load modules
    we are just interested in the first domains on a gene.
    when we encounter a new gene,
    we want to stop checking if the gene has a starter module.
    we also want to start checking again when we go over to a new gene.

    '''
    # First we check to see if antismash has located a typical A/AT-ACP/PCP-module:
    # if this is the case, we dont want to add a loader module.

    for module in domain_or_module:
        if _is_standard_load_module(module):
            return domain_or_module

    found_loader = False
    prev_gene = ''
    encountered_a_module = False
    for index, dom_mod in enumerate(domain_or_module):  # dom_mod is either a domain or a module
        if dom_mod['element'] == 'domain':  # if this is a domain and not a module
            if dom_mod['info']['gene'] != prev_gene:
                encountered_a_module = False
                start_of_gene = index  
                # so we know what domains to remove in case we find a loader domain
                # if we are on a new gene, the search continues. 
                # this is only relevant to do on domains, as modules cannot (should not) by antismash definition
                # should always be extending, and should in theory not contain cal,fkbh or gnat domains
            if dom_mod['info']['core_gene'] and not encountered_a_module:
                # if this is a core gene. needs to be nested because modules have different structure than domains here
                # Check if the domain is either Cal, GNAT or FkbH
                if dom_mod['info']['type'] in loader_domains:  
                    # if this is a gene that is a giveaway that this is the loader module:
                    # we have found a loader domain on one of the core genes, before any other modules have been found.
                    found_loader = True
                    loader_type = dom_mod['info']['type']
                    loader_activity = dom_mod['info']['activity']
                    loader_start = start_of_gene
                    loader_gene = dom_mod['info']['gene']

            prev_gene = dom_mod['info']['gene']
            if found_loader and dom_mod['info']['type'] in loader_acp_domains:
                new_domains = [x['info'] for x in domain_or_module[loader_start:index+1]]
                info = {'extender_unit': loader_type,
                        'loader_activity': loader_activity,
                        'start': domain_or_module[start_of_gene]['info']['start'],
                        'end': domain_or_module[index]['info']['end']}
                new_module = {'element': 'module', 'info': info,'domains': new_domains}
                return domain_or_module[:start_of_gene] + [new_module] + domain_or_module[index+1:]

        elif dom_mod['element'] == 'module':
            # to check that we have not passed any real modules
            encountered_a_module = True  
            # If we have previously found either a Cal, GNAT or FkbH loader domain
            # and this module is on the same gene we have now found the starter module
            if found_loader and dom_mod['domains'][0]['gene'] == prev_gene:
                info = {'extender_unit': loader_type,
                        'loader_activity': loader_activity,
                        'start': domain_or_module[start_of_gene]['info']['start'],
                        'end': domain_or_module[index]['info']['end']}
                new_domains = [x['info'] for x in domain_or_module[loader_start:index-1]]
                new_module = {'element': 'module', 'info': info,'domains': new_domains}
                return domain_or_module[0:start_of_gene] + [new_module] + domain_or_module[index:]

            if not dom_mod['domains'][0]['gene'] == prev_gene:
                # so we know what domains to remove in case we find a loader domain
                start_of_gene = index  

            # the [0] part is just because modules can in theory span across several genes
            prev_gene = dom_mod['domains'][0]['gene']
            

    # Then we check if there are any loader modules that can be found within the core genes that are not inside a
    # module as defined by antiSMASH. These are the [A/AT]-[ACP/PCP/PP] modules.
    prev_domain_type = ''
    prev_act = ''
    prev_gene = ''
    # encountered_a_module = False
    for index, dom_mod in enumerate(domain_or_module):  # dom_mod is either a domain or a module
        if dom_mod['element'] == 'domain':  # if this is a domain and not a module (implicitly not a module)
            if not dom_mod['info']['gene'] == prev_gene:  # If the domain is the first domain on a gene.
                started_new_gene = True

            # starter_module
            if dom_mod['info']['core_gene']:# and not encountered_a_module:
                # if this is a core gene. needs to be nested because modules have different structure than domains
                if prev_domain_type in loader_at_domains and dom_mod['info']['type'] in loader_acp_domains:
                    if prev_gene == dom_mod['info']['gene']:
                        # if this is a gene that contains the [A/AT]-[ACP/PCP/PP] module outside of a regular
                        # extending module:
                        # try:
                        # print(domain_or_module[index - 1])
                        extender_unit = domain_or_module[index - 1]["info"]["AT_specificity"]
                        # except KeyError:
                        #     extender_unit = 'starter_unit'

                        info = {'extender_unit': extender_unit,
                                'activity': prev_act,
                                'start': domain_or_module[index-1]['info']['start'],
                                'end': domain_or_module[index]['info']['end']}
                        new_domains = [domain_or_module[index - 1]['info'], domain_or_module[index]['info']]
                        new_module = {'element': 'module', 'info': info,'domains': new_domains}
                        return domain_or_module[:index] + [new_module] + domain_or_module[index:]

            prev_domain_type = dom_mod['info']['type']
            prev_act = dom_mod['info']['activity']
            prev_gene = dom_mod['info']['gene']


        elif dom_mod['element'] == 'module':
            # Only looking for loader domains before the modules
            break
            # encountered_a_module = True  # to check that we have not passed any real modules
            
            # prev_gene = dom_mod['domains'][0]['gene']

    # If we still haven't found a load module we check if the load module is a C-A-PCP module for NRPS
    prev_domain_types = ["", "", ""]
    for index, dom_mod in enumerate(domain_or_module):
        # Update the prev domain list
        if dom_mod['element'] == 'domain':
            prev_domain_types.pop(0)
            prev_domain_types.append(dom_mod['info']['type'])

            if prev_domain_types == NRPS_acylating_loader:
                # We have found a particular C-A-PCP loader module
                info = {'extender_unit': 'NRPS_acylating_loader',
                        'activity': prev_act,
                        'start': domain_or_module[index]['info']['start'],
                        'end': domain_or_module[index]['info']['end']}
                new_domains = [domain_or_module[index-2]]
                new_module = {'element': 'module', 'info': info,'domains': new_domains}

                return domain_or_module[:index] + [new_module] + domain_or_module[index:]

        elif dom_mod['element'] == 'module':
            # Either the first module has these three domains or we don't know what the load module is
            if _all_NRPS_acylating_domains(dom_mod):
                
                # Have found a loader C-A-PCP module
                info = {'extender_unit': 'NRPS_acylating_loader',
                        'activity': prev_act,
                        'start': domain_or_module[index]['info']['start'],
                        'end': domain_or_module[index]['info']['end']}
                new_domains = [dom_mod["domains"][0]]
                new_module = {'element': 'module', 'info': info,'domains': new_domains}
                return domain_or_module[:index] + [new_module] + domain_or_module[index:]
            else:
                break


    # in case we dont find any loader module, we assume its existence anyways, based on the finding that
    # this is usually the case. We first try to find a free-standing A or AT-domain in order to predict wether the
    # Starter is an amino acid or a Acyl-unit.
    domain_or_module.insert(0, {'element': 'module', 'info': {
        'extender_unit': 'custom_starter', 'start': 0,
        'end': 0}, 'domains': []})

    return domain_or_module

def _all_NRPS_acylating_domains(module):
    module_domain_types = [d["type"] for d in module["domains"]]
    for d in NRPS_acylating_loader:
        if not d in module_domain_types:
            return False
    return True


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

    
    domains, actual_domains, domains_x_modules = _get_domains(data, core_number)

    # assume that DHD-domains are inactive, however, may still dehydratase the previous domain(?)
    # assue that all others are active.

    domains_x_modules = find_and_replace_load_modules(domains_x_modules)
    domains_x_modules = deactivate_reducing_domains_when_kr_is_inactive(domains_x_modules, True)
    model = create_t1_transat_nrps_model(data['core_structure'][core_number], domains_x_modules, model,
                                         tailoring_reactions)


    # domains x modules is now a dict sorted by placement on BGC.
    # do this to find out if domains are found in this sequence: ks-at-acp
    return model


def _get_domains(data, core_number):
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

    return domains, actual_domains, domains_x_modules

def add_nrps_metabolic_pathway(data, model, core_number, tailoring_reactions):
    """
    :param data: all data that has been structured from gbk-file
    :return: modules max: The maximum detected modules, following rule that it needs a KS-domain and an ACP domain
    Modules mid: All modules that contain KS, AT and ACP domain.
    the modules_min is
    """
    
    domains, actual_domains, domains_x_modules = _get_domains(data, core_number)

    domains_x_modules = find_and_replace_load_modules(domains_x_modules)
    model = create_t1_transat_nrps_model(data['core_structure'][core_number], domains_x_modules, model,
                                         tailoring_reactions)


    return model


def force_X_nrps_module_flux(substrate):
    '''
    :param substrate: an amino acid that is either not specifically determined (NRPS predictor couldnt find
    specific amino acid)
    :return: 20 reactions that convert each amino acid into the unknown amino acid
    '''
    met_id = long_to_bigg[substrate]
    metabolite = _new_met(met_id, name='Generic amino acid')
    reaction_list_aa = []
    for amino_acid in acid_to_bigg:
        reaction = cobra.Reaction(amino_acid + '_to_' + substrate)
        reaction.name = 'Convert to generic metabolite'
        reaction.lower_bound = 0  # This is the default
        reaction.upper_bound = 1000
        reaction.add_metabolites({ref_model.metabolites.get_by_id(acid_to_bigg[amino_acid]): -1, metabolite: 1})
        reaction_list_aa.append(reaction)
    return reaction_list_aa, metabolite

def _add_metabolites(r, met_dict):
    for key, value in met_dict.items():
        if isinstance(key, str):
            m = cobra.Metabolite(key, formula = _get_aa_formula(key))
        else:
            m = key
        r.add_metabolites({m: value})


def _make_methoxymalCoA_reaction():
    mx_rx = cobra.Reaction('mxmal_synthesis')
    mx_rx.name = 'synthesis of methoxymalonyl-coa'
    mx_rx.lower_bound = 0.  # This is the default
    mx_rx.upper_bound = 1000.
    mx_rx.add_metabolites(cofactor_reactions_dict["mxmal"])
    return mx_rx

def _make_IN_reaction(core_structure):
    in_rx = cobra.Reaction('IN_' + core_structure['type'])
    in_rx.name = core_structure['type'] + '_reaction'
    in_rx.lower_bound = 0.  # This is the default
    in_rx.upper_bound = 1000.
    new_met = _new_met(core_structure['type'] + '_0', name=core_structure['type'])
    in_rx.add_metabolites({new_met: 1.0})
    return in_rx

def _new_met(met_id, name = None, formula = "X"):
    if not name:
        name = met_id
    new_met = cobra.Metabolite(met_id, formula = formula, name  = name, compartment = "c")
    return new_met



def _get_gene_reaction_rule(module):
    genes = []
    for d in module["domains"]:
        try:
            gene = d["gene"][0]
        except KeyError:
            continue
        else:
            if len(gene):
                genes.append(gene)

    genes = list(set(genes))
    if len(genes) == 0:
        return ""
    elif len(genes) == 1:
        return genes[0]
    else:
        return " and ".join(genes)

def _get_module_type(module):
    """
    Determine the type of module (PKS or NRPS) based on the domains in the module
    """
    NRPS_signature_domains = ["Condensation", "AMP-binding"]
    PKS_signature_domains = ["PKS_KS", "PKS_AT"]
    if module['info']['extender_unit'] == 'custom_starter':
        return 'PKS'
    else:
        for domain in module['domains']:
            if domain['type'] in PKS_signature_domains:
                return 'PKS'

            elif domain['type'] in NRPS_signature_domains:
               return 'NRPS'
            else:
                pass

    # No signature domains. Look for specific starters
    for domain in module['domains']:
        if domain['type'] in alternate_starters:
            return 'starter_module'

    return None

def create_t1_transat_nrps_model(core_structure, domains_x_modules, model, tailoring_reactions):
    # vi lager modellen som inneholder alle reaksjonene.
    reaction_list = []  # list of reactions that take place
    if tailoring_reactions['methoxymalonyl']:
        reaction_list.append(_make_methoxymalCoA_reaction())

    # Make an sink reaction for the initial start of the polyketide, this is required in the following loop
    in_rx = _make_IN_reaction(core_structure)

    # Flag to determine if the polymer is cleaved off by a TE domain
    chain_released = False
    
    # Make a flag so the conversion reactions from specific to generic AA's only are added once
    domain_counter = 1
    for module in domains_x_modules:
        # We are only interested in modules
        if module['element'] == 'domain':
            continue
        
        module_type = _get_module_type(module)

        if module_type in ["PKS", "NRPS", "starter_module"]:
            # This needs to be first because we want to add extender unit as first step of a module, but for transAT
            # modules, the at domain is not always present.
            if (module['info']['extender_unit'] not in non_extending_modules) and not chain_released:
                # essentially: if the module is a regular extender module
                # Get extender unit
                # if module_type == 'starter_module':
                #     extender_unit = return_key_by_value(module['info']['activity'].split(' -> ')[0])
                # else:
                extender_unit = return_key_by_value(module['info']['extender_unit'].split(' -> ')[0])

                # This is just a way that we can separate extender units that exist in the GEM
                # from those that do not exist in the GEM. if the TRY fails, it means that
                # the extender unit does not exist in the GEM
                # happens if extender unit is not malcoa, mmcoa, mxcoa or emcoa,
                # (or minowa and at_specificity disagrees):
                if extender_unit == 'pk':

                    # if this is the empty loader we created:
                    if module['info']['extender_unit'] == 'custom_starter':
                        extender_unit = 'Malonyl-CoA'

                    # if this is an alternate loader, i.e gnat, fkbh, cal loader:
                    elif module['info']['extender_unit'] in alternate_starters:
                        if module['info']['extender_unit'] == 'CAL_domain':
                            module_type = module['info']['loader_activity']
                            print("#!#", module['info']['loader_activity'])
                        elif module['info']['extender_unit'] == 'FkbH':
                            module_type = 'FkbH'
                        elif module['info']['extender_unit'] == 'GNAT':
                            module_type = 'GNAT'
                        elif module['info']['extender_unit'] == 'starter_unit':
                            module_type = 'starter_module'
                        elif module['info']['extender_unit'] == 'NRPS_acylating_loader':
                            module_type = 'NRPS_acylating_loader'


                    # then  we check if this should be mxmal:
                    # if minowa and at_prediction had no consensus, we set extender unit to mxmal if there is
                    # machinery that synthesises it. if not, we assume the extender is malcoa
                    for domains in module['domains']:
                        if domains['type'] == 'PKS_AT':
                            if tailoring_reactions['methoxymalonyl'] and domains['minowa'] == 'mxmal' or domains[
                                'AT_specificity'] == 'mxmal':
                                extender_unit = 'Methoxymalonyl-CoA'
                            else:
                                # we have to choose between minowa and AT_specificity, so we chose minowa:
                                # this also has the consequence that if minowa is emal: the next step may convert
                                # the extender into methoxymalonyl-CoA
                                extender_unit = return_key_by_value(domains['minowa'])

                if tailoring_reactions['methoxymalonyl'] and extender_unit == 'Ethylmalonyl-CoA':
                    extender_unit = 'Methoxymalonyl-CoA'

                prevmet = _new_met(core_structure['type'] + '_' + str(domain_counter - 1), core_structure['type'])
                postmet = _new_met(core_structure['type'] + '_' + str(domain_counter), core_structure['type'])
                
                reaction = cobra.Reaction(core_structure['type'] + '_' + str(domain_counter))
                reaction.name = core_structure['type'] + '_reaction_' + str(domain_counter)
                reaction.lower_bound = 0.  # This is the default
                reaction.upper_bound = 1000.
                reaction.gene_reaction_rule = _get_gene_reaction_rule(module)

                if module_type == 'NRPS':
                    if extender_unit == "hpg":
                        hpg_synthesis_reactions = make_hpg_reactions()
                        model.add_reactions(hpg_synthesis_reactions)
                        extender_met = cofactor_metabolites_dict["hpg"]
                    elif extender_unit == 'dpg' or extender_unit == 'dhpg':
                        # duplicate entries
                        # dpg and dhpg are the same substrates
                        dpg_reaction1 = cobra.Reaction('dhpg_synthesis')
                        dpg_reaction1.name = 'synthesis of dhpg'
                        dpg_reaction1.lower_bound = 0.  # This is the default
                        dpg_reaction1.upper_bound = 1000.
                        dpg_reaction1.add_metabolites(cofactor_reactions_dict['dpg'])
                        model.add_reactions([dpg_reaction1])
                        extender_met = cofactor_metabolites_dict["dpg"]
                        
                    elif extender_unit in ['pip', 'bht', 'abu']:
                        synthesis_reaction = cobra.Reaction('{0}_synthesis'.format(extender_unit))
                        synthesis_reaction.name = 'synthesis of {0}'.format(std_aa_dic[extender_unit])
                        synthesis_reaction.lower_bound = 0.  # This is the default
                        synthesis_reaction.upper_bound = 1000.
                        synthesis_reaction.add_metabolites(cofactor_reactions_dict[extender_unit])
                        model.add_reactions([synthesis_reaction])
                        extender_met = cofactor_metabolites_dict[extender_unit]

                    else:
                        try:
                            extender_met = ref_model.metabolites.get_by_id(long_to_bigg[extender_unit])
                        except KeyError:
                            # There are a range of very specific metabolites that are not accounted for
                            # including the cases where antiSMASH can't predict the specific AA
       
                            AA_to_X_reactions, extender_met = force_X_nrps_module_flux(extender_unit)
                            model.add_reactions(AA_to_X_reactions)

                    
                    reaction.add_metabolites(
                        {extender_met:-1,
                         ref_model.metabolites.get_by_id('atp_c'): -1,
                         ref_model.metabolites.get_by_id('amp_c'): 1,
                         ref_model.metabolites.get_by_id('ppi_c'): 1,
                         ref_model.metabolites.get_by_id('h2o_c'): 1})
                    
                elif module_type == 'PKS':
                    if extender_unit == 'Methoxymalonyl-CoA' or extender_unit == 'Ethylmalonyl-CoA':
                        extender_met = _new_met(long_to_bigg[extender_unit], name = extender_unit)
                    else:
                        try:
                            extender_met = ref_model.metabolites.get_by_id(long_to_bigg[extender_unit])
                        except KeyError:
                            print(extender_unit, " is not in model, assume malonyl-CoA")
                            # If we don't now what the extender unit is we assume malonyl-CoA
                            extender_met = ref_model.metabolites.get_by_id(long_to_bigg["Malonyl-CoA"])

                    reaction.add_metabolites({extender_met: -1, ref_model.metabolites.get_by_id('coa_c'): 1,
                                              ref_model.metabolites.get_by_id('co2_c'): 1})
                    
                elif module_type == 'FkbH':
                    reaction.add_metabolites(cofactor_reactions_dict['fkbh'])
                    
                elif module_type == 'GNAT':
                    reaction.add_metabolites(cofactor_reactions_dict['gnat'])
                    
                elif module_type == 'fatty_acid':
                    reaction.add_metabolites(cofactor_reactions_dict['fatty_acid'])
                    
                elif module_type == 'AHBA':
                    ahba_synthesis_reaction = cobra.Reaction('ahba_synthesis')
                    ahba_synthesis_reaction.name = "AHBA synthesis"
                    ahba_synthesis_reaction.bounds = (0,1000)
                    ahba_synthesis_reaction.add_metabolites(cofactor_reactions_dict['ahba-synthesis'])
                    model.add_reaction(ahba_synthesis_reaction)
                    reaction.add_metabolites(cofactor_reactions_dict['ahba'])
                    
                elif module_type == 'shikimic_acid':
                    reaction.add_metabolites(cofactor_reactions_dict['shikimic_acid'])
                    
                elif module_type == 'Acetyl-CoA':
                    reaction.add_metabolites(cofactor_reactions_dict['acetyl'])
                    
                elif module_type == 'NH2':
                    reaction.add_metabolites(cofactor_reactions_dict['NH2'])
                elif module_type == "NRPS_acylating_loader":
                    reaction.add_metabolites(cofactor_reactions_dict["NRPS_acylating_loader"])
                    # Add reactions making a generic fatty acid
                    fatty_acid_generic_reactions = make_generic_fatty_acid_conversion()
                    model.add_reactions(fatty_acid_generic_reactions)
                    
                elif module_type == 'starter_module':
                    print("Module type is starter module")
                    raise ValueError
                    # reaction.add_metabolites(
                    #     {prevmet: -1,
                    #      postmet: 1})
                    # reaction.add_metabolites(cofactor_reactions_dict['NH2'])
                    # reaction_list.append(reaction)
                    # domain_counter += 1
                else:
                    # raise arbitrary error because something has gone terribly wrong
                    raise ValueError

                # Add the polyketide chain metabolites to the reation and add reaction to list
                
                # remove CO2 from load module
                if (domain_counter == 1) and (module_type != 'GNAT'):
                    remove_one_co2(reaction)

                reaction.add_metabolites({prevmet: -1, postmet: 1})
                reaction_list.append(reaction)
                domain_counter += 1


                # except KeyError:  # the substrate is not in the model you want to insert into
                    
                    # warnings.warn('Predicted substrate does not exist in the model')
                    # extender_unit = return_key_by_value(module['info']['extender_unit'].split(' -> ')[0])

                    # if module_type == 'PKS':
                    #     extender_unit = 'Malonyl-CoA'
                        # the only cases we end up here when the module type is PKS, is when a
                        # starter unit is predicted as the extender unit. This Cannot be true for starters that antiSMASH
                        # predicts, so we set this extender to malonyl-CoA

                    # prevmet = _new_met(core_structure['type'] + '_' + str(domain_counter - 1), core_structure['type'])
                    # postmet = _new_met(core_structure['type'] + '_' + str(domain_counter), core_structure['type'])

                    # # Make reaction
                    # reaction = cobra.Reaction(core_structure['type'] + '_' + str(domain_counter)) 
                    # reaction.name = core_structure['type'] + '_reaction_' + str(domain_counter)
                    # reaction.lower_bound = 0.  # This is the default
                    # reaction.upper_bound = 1000.
                    # if module_type == 'NRPS':
                        

                        
                    #     else:
                    #         print(module_type)
                    #         print(module)
                    #         # we dont know the specificity of the A-domain, so we let any amino acid take the role
                    #         # as this amino acid. 
                    #         if not added_generic_AAs:
                    #             model.add_reactions(force_X_nrps_module_flux(extender_unit))
                    #             added_generic_AAs = True
                    #         new_met = cobra.Metabolite(long_to_bigg[extender_unit], formula='X', name=extender_unit,compartment='c')
                    #         reaction.add_metabolites({
                    #             new_met: -1,
                    #             ref_model.metabolites.get_by_id('atp_c'): -1,
                    #             ref_model.metabolites.get_by_id('amp_c'): 1,
                    #             ref_model.metabolites.get_by_id('ppi_c'): 1,
                    #             prevmet: -1,
                    #             postmet: 1})
                    #     reaction_list.append(reaction)
                    #     warnings.warn('Predicted substrate does not exist in the model')
                    #     domain_counter += 1
                    # if module_type == 'PKS':
                    #     reaction.add_metabolites(
                    #         {cobra.Metabolite(long_to_bigg[extender_unit], formula='X', name=extender_unit,
                    #                           compartment='c'): -1,
                    #          ref_model.metabolites.get_by_id('coa_c'): 1,
                    #          ref_model.metabolites.get_by_id('co2_c'): 1,
                    #          prevmet: -1,
                    #          postmet: 1})
                    #     reaction_list.append(reaction)
                    #     domain_counter += 1
                    # else:
                    #     raise IndexError
            
            for domain in module['domains']:
                if domain['activity']:
                    if domain['type'] in general_domain_dict:
                        std_domain_name = general_domain_dict[domain['type']]
                        reaction = cobra.Reaction(core_structure['type'] + '_' + str(domain_counter))
                        reaction.name = std_domain_name
                        reaction.lower_bound = 0.  # This is the default
                        reaction.upper_bound = 1000.
                        reaction.add_metabolites(cofactor_reactions_dict[std_domain_name])

                        prevmet = _new_met(core_structure['type']+'_'+ str(domain_counter - 1), core_structure['type'])
                        postmet = _new_met(core_structure['type']+'_'+ str(domain_counter), core_structure['type'])
                        reaction.add_metabolites({prevmet: -1, postmet: 1})
                        reaction_list.append(reaction)
                        domain_counter += 1

                        # If we encounter a TE domain the polymer is released and it is unlikely that 
                        # further elongation occurs
                        if domain["type"] in chain_release_domains:
                            chain_released = True

        # else:
        #     print(module)
        #     raise ValueError

    '''
    this part below is to add tailoring reactions:
    tailoring reactions is a dict of type:

    {glycosyltransferase: 3, ALA: 0, glycerol: 1}
    '''
    for tailoring_reaction, n in tailoring_reactions.items():
        if tailoring_reaction != 'methoxymalonyl':
            for repetition in range(n):
                reaction = cobra.Reaction(core_structure['type'] + '_' + str(domain_counter))
                reaction.add_metabolites(tailoring_reactions_dict[tailoring_reaction])
                prevmet = _new_met(core_structure['type'] + '_' + str(domain_counter - 1), core_structure['type'])
                postmet = _new_met(core_structure['type'] + '_' + str(domain_counter), core_structure['type'])
                reaction.add_metabolites({prevmet: -1.0, postmet: 1.0})
                reaction_list.append(reaction)
                domain_counter += 1

    '''
    The reaction below is here in order to have a reaction that converts e.g. T1PKS_54 to a general metabolite
    we do this so that we can add tailoring reactions. 
    '''
    ex_rx = cobra.Reaction('DM_secondary_metabolite')
    ex_rx.name = core_structure['type'] + '_reaction'
    ex_rx.lower_bound = 0.  # This is the default
    ex_rx.upper_bound = 1000.
    ex_rx.add_metabolites({postmet: - 1.0})

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

def remove_one_co2(reaction):
    """
    CO2 is only produced in PKS chain elongations (claissen condensations). Thus for the load module CO2 is not produced
    """
    for m, i in reaction.metabolites.items():
        if m.id == "co2_c":
            if i > 0:
                reaction.add_metabolites({m:-1})
                break
        

def add_ripp_metabolic_pathway(core_structure, model, core_number):
    reaction_name = "{0}_{1}".format(core_structure['type'], core_number)
    reaction = cobra.Reaction(reaction_name)
    reaction.name = core_structure['type'] + '_reaction'
    reaction.lower_bound = 0.  # This is the default
    reaction.upper_bound = 1000.

    ripp_met = _new_met(reaction_name, name=core_structure['type'])
    reaction_metabolites = {ripp_met: 1.0}
    aa_metabolites = {}
    for letter in core_structure['RiPP']:
        if letter in aa_metabolites:
            aa_metabolites[letter] -= 1
        else:
            aa_metabolites[letter] = -1

    for metabolite in aa_metabolites:  # color red and metabolites_snorre

        try:
            met_bigg = long_to_bigg[one_letter_aa[metabolite]]
            metabolite_added = ref_model.metabolites.get_by_id(met_bigg)
            reaction_metabolites[metabolite_added] = aa_metabolites[metabolite]
        except KeyError:
            continue

    reaction.add_metabolites(reaction_metabolites)

    # Make a generic end product (to handle multiple products)
    generic_met_id = "generic_{0}".format(core_structure['type'])
    try:
        generic_met = model.metabolites.get_by_id(generic_met_id)
    except KeyError:
        generic_met = _new_met(generic_met_id, name=core_structure['type'])

    translation_reaction_name = "translate_{0}".format(core_number)
    translate_reaction = cobra.Reaction(translation_reaction_name)
    translate_reaction.add_metabolites({ripp_met:-1, generic_met:1})

    model.add_reactions([reaction, translate_reaction])

    # Add DM reaction for RiPP
    dm_reaction_name = "DM_secondary_metabolite"
    try:
        r = model.reactions.get_by_id(dm_reaction_name)
    except KeyError:
        ex_rx = cobra.Reaction(dm_reaction_name)
        ex_rx.name = core_structure['type'] + '_reaction'
        ex_rx.lower_bound = 0.  # This is the default
        ex_rx.upper_bound = 1000.
        ex_rx.add_metabolites({generic_met: - 1.0})
        model.add_reactions([ex_rx])

    return model

def make_generic_fatty_acid_conversion():
    reactions = []
    for name, met in fatty_acyl_CoAs.items():
        reaction = cobra.Reaction("{0}_to_generic".format(name))
        reaction.lower_bound = 0
        reaction.upper_bound = 1000
        mets = {met:-1, cofactor_metabolites_dict["fatty_acid_X"]:1}
        reaction.add_metabolites(mets)
        reactions.append(reaction)
    return reactions

def make_hpg_reactions():
    """
    Make reactions required to syntesize the extender unit hpg
    """
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

    return [hpg_reaction1, hpg_reaction2, hpg_reaction3, hpg_reaction4]

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

    # FkbH_like protein with smcog number 1256 - for methoxymalonyl-coa
    tailoring_genes_based_on_smCOG_definition_dict[1256] = False  
    # 3-Hydroxybutyryl-CoA dehydrogenase like protein with smcog number 1095 - for methoxymalonyl-coa
    tailoring_genes_based_on_smCOG_definition_dict[1095] = False  
    # glycosyltransferase with smcog 1062 - For glycosyl groups
    tailoring_genes_based_on_smCOG_definition_dict[1062] = False  
    # for tailoring reaction that adds glycerol to PK through 1,3-biphosphoglycerate
    tailoring_genes_based_on_smCOG_definition_dict[1084] = False  
    # AMP-dependent synthase and ligase with smcog 1002 - 2-Amino-3-hydroxycyclopent-2-enone
    tailoring_genes_based_on_smCOG_definition_dict[1002] = False  
    # 8-amino-7-oxononanoate synthase with smcog 1109- For 2-Amino-3-hydroxycyclopent-2-enone
    tailoring_genes_based_on_smCOG_definition_dict[1109] = False  

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


def add_cores_to_model(name, data, model_output_path):
    model = cobra.Model(name)

    smcog_dict = find_tailoring_reactions_from_smcogs(data)

    tailoring_reactions_dict = add_tailoring_smcogs(smcog_dict)

    bgc_types = []
    for core_number in data['core_structure']:
        print("Core number: ", core_number)
        bgc_type = data['core_structure'][core_number]['type']
        bgc_types.append(bgc_type)
        
        if bgc_type in RiPPs:
            print("Ripps not successfully implemented yet")
            pass
            #model = add_ripp_metabolic_pathway(data['core_structure'][core_number], model, core_number)

        elif bgc_type in ['transAT-PKS', 'transAT-PKS-like']:
            model, lump_model = add_transat_metabolic_pathway(data, model, core_number, tailoring_reactions_dict)

        elif bgc_type == 'T1PKS':
            model, lump_model = add_t1pks_metabolic_pathway(data, model, core_number, tailoring_reactions_dict)

        elif bgc_type in ['NRPS', 'NRPS-like']:
            model, lump_model = add_nrps_metabolic_pathway(data, model, core_number, tailoring_reactions_dict)
        else:
            print("Can't handle ", bgc_type)
    '''
    Then finally, we save the model to a file
    '''
    if not len(model.reactions):
        print("Couldn't construct the pathway for ", name)
        print(smcog_dict)
        print(data['core_structure'])
        return False, "-".join(bgc_types)
    else:
        model.description =  "-".join(bgc_types)
        cobra.io.save_json_model(model, model_output_path)
        return True, "-".join(bgc_types)

def run(bgc_path, output_folder, json_folder = None):
    """
    The main function used to run BiGMeC on wither a single .gbk file or a folder
    """
    bgc_path = Path(bgc_path)
    report_list = []
    if bgc_path.is_dir():
        for filename in bgc_path.glob("*.gbk"):
            print(filename)
            #Run this script for each file in the folder
            result = _run(filename, output_folder, json_folder)
            report_list.append(result)
    else:
        result = _run(bgc_path, output_folder, json_folder)
        report_list.append(result)

    df = pd.DataFrame(report_list, columns = ["BGC", "Success", "BGC type"])
    df.to_csv(output_folder + "/summary.csv")

def _run(bgc_path, output_folder, json_folder):
    # This is the core of this function
    json_path = str(Path(json_folder) / (bgc_path.stem + '.json'))
    data = parse_antismash_gbk(bgc_path, json_path)
    
    # Adds extracted data to model
    output_model_fn = Path(output_folder) / (bgc_path.stem + ".json")
    name = "BGC-{0}".format(bgc_path.stem)
    successfull, _ = add_cores_to_model(name, data, str(output_model_fn))
    bgc_type_str = "/".join(data["Cluster types"])
    return bgc_path.stem, int(successfull), bgc_type_str

if __name__ == '__main__':
    biggbk = "../Data/mibig"
    
    # 1) Folder containing all gbk files you want to translate into metabolic pathways
    #    They are saved as models, and can be merged with the GEM later.
    #    In this repository, the models that are found in "gbk_db_output_models.zip" are the pathways that 
    #    Have been constructed. In that case, this folder contained all antiSMASH output files found in 
    #    antismash_output_mibig_gbk_files1.zip
    #    antismash_output_mibig_gbk_files2.zip
    #    and
    #    antismash_output_mibig_gbk_files3.zip
    output_gbk = "../Data/constructed_pathways/"
    # 2) Folder that regular models are output (empty folder)
    output_gbk_lump = "../Data/constructed_GEMs"
    # 3) Folder that lump models are output (empty folder)
    json_folder = '../Data/temp'
    # 4) Folder that json files are output (empty folder)(an unnecessary step, but it may help look at the information that is saved)

    # 5) path of the genome scale metabolic model
    # model_fn = '../Models/Sco-GEM.xml'
    # ref_model = cobra.io.read_sbml_model(model_fn)

    if 0:
        for filename in os.listdir(biggbk):
            if filename.endswith('.gbk'):
                json_path = json_folder + filename.split('.')[0] + '.json'
                parse_antismash_gbk(biggbk + filename)
                with open(json_path, 'r') as json_file:
                    data_json = json.load(json_file)
                json_file.close()
                add_cores_to_model(data_json, output_gbk + filename[:-4] + ".json")
    if 1:

        bgc_path = biggbk #+ "/1049.gbk"
        run(bgc_path, output_gbk, json_folder)
