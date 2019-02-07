#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This is untested

from scripts.ot_genetics_client import OT_Genetics

# Get client
otg = OT_Genetics()

# List of rsids
rsid_list = ['rs123', 'rs1234']

# Convert rsids to variant IDs.
# Returns dictionary, e.g. {'rs123': ['1_3423_A_T', '1_3423_A_G']}
rsid_to_varid_dict = otg.rsid_to_varid(rsid_list)

# Flatten variant ID lists
def flatten(l):
    return [item for sublist in l for item in sublist]
varid_list = flatten(rsid_to_varid_dict.values())

# Get V2G data. Returns pandas dataframe
v2g = otg.get_v2g_evidence_for_varids(varid_list)
