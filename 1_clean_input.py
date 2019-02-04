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
import numpy as np

def main():

    # Args
    in_pattern = 'validation_sets/raw/*.tsv'
    out_dir = 'validation_sets/cleaned'

    os.makedirs(out_dir, exist_ok=True)

    # Get client
    otg = OT_Genetics()

    # Process all input sets
    for inf in glob(in_pattern):

        # Check if output exists
        outf = os.path.join(out_dir, os.path.split(inf)[1])
        if os.path.exists(outf):
            print('Warning, {} exists, skipping...'.format(outf))
            continue

        # Clean the input
        clean = process_raw(inf, otg)

        # Write
        clean.to_csv(outf, sep='\t', index=None)


    return 0

def process_raw(inf, client):
    ''' Loads and cleans the raw into validation set
    Params:
        inf (str): input file
        client (OT_Genetics)
    Returns:
        df
    '''
    # Load
    df = pd.read_csv(inf, sep='\t', header=0)

    # Make assertations
    assert (df.rsid.str.startswith('rs')).all()
    assert (df.condifence.isin(['High', 'Medium', 'Low'])).all()
    assert ~(pd.isnull(df.rsid) & pd.isnull(df.varid)).all()
    assert ~pd.isnull(df.loc[:, ['set', 'source', 'gene_id', 'condifence']]).any(axis=None)

    # Find unkown variant IDs
    is_unknown = pd.isnull(df.varid)
    rsids_to_get = df.rsid[is_unknown]
    rsid_map = client.rsid_to_varid(rsids_to_get.tolist())
    df.loc[is_unknown, 'varid'] = df.loc[is_unknown, 'rsid'].apply(
        lambda rsid: rsid_map.get(rsid, [])
    )
    df['varid_ambiguous'] = (df.varid.str.len() > 1)

    # Explode lines
    df = explode(df, ['varid'])

    return df

def explode(df, columns):
    ''' Explodes multiple columns
    '''
    idx = np.repeat(df.index, df[columns[0]].str.len())
    a = df.T.reindex(columns).values
    concat = np.concatenate([np.concatenate(a[i]) for i in range(a.shape[0])])
    p = pd.DataFrame(concat.reshape(a.shape[0], -1).T, idx, columns)
    return pd.concat([df.drop(columns, axis=1), p], axis=1).reset_index(drop=True)

if __name__ == '__main__':

    main()
