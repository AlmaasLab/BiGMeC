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