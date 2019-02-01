#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
from glob import glob
import pandas as pd
from scripts.ot_genetics_client import OT_Genetics

def main():

    # Args
    in_pattern = 'validation_sets/cleaned/*.tsv'
    out_dir = 'data/v2g/190201'

    os.makedirs(out_dir, exist_ok=True)

    # Get client
    otg = OT_Genetics()

    # Get set of variant IDs
    varids = set([])
    for inf in glob(in_pattern):
        df = pd.read_csv(inf, sep='\t', header=0)
        varids = varids.union(df.varid)

    # Query API
    res = otg.get_v2g_evidence_for_varids(list(varids))

    # Write
    outf = os.path.join(out_dir, 'evidence.tsv')
    res.to_csv(outf, sep='\t', index=None)

    return 0

if __name__ == '__main__':

    main()
