#!/usr/bin/env python
# coding: utf-8
"""
Copyright 2020 Snorre Sulheim (snorre.sulheim@sintef.no)
https://github.com/AlmaasLab/BiGMeC

This file holds lists of domain names that are used in the bigmec pipeline


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
NRPS_acylating_loader = ["Condensation", "AMP-binding", "PCP"]

non_extending_modules = [
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
    "nMT": "nMT",
    "Thioesterase": "PKS_TE",
    "TD": "TD"
}

acp_domains = ["ACP", "PP-binding", 'ACP_beta']
at_domains = ['Trans-AT_docking', 'PKS_AT']
trans_at_cores = ['transAT-PKS-like', 'transAT-PKS']
nrps_pks_cores = ['T1PKS', 'NRPS', 'PKS-like', 'NRPS-like']
compatible_cores = ['transAT-PKS-like', 'transAT-PKS', 'T1PKS', 'NRPS', 'PKS-like', 'NRPS-like']
RiPPs = ['lanthipeptide', 'thiopeptide', 'lassopeptide']
prefer_types = ['transAT-PKS', 'transAT-PKS-like', 'NRPS', 'T1PKS', 'NRPS-like', 'PKS-like']
dh_domains = ['PKS_DH', 'PKS_DH2', "PKS_DHt",]
print_list = ['PKS_KS', 'PKS_DH', 'PKS_DH2', 'PKS_KR', 'MT', 'PKS_ER', 'Condensation', 'Thioesterase']
alternate_starters = ['CAL_domain', 'FkbH', 'GNAT', 'NRPS_acylating_loader']
chain_release_domains = ["Thioesterase", "TD"]