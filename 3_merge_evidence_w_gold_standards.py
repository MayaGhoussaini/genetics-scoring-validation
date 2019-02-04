#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import pandas as pd

def main():

    # Args
    in_gold_standards = 'validation_sets/cleaned/progem_190201.tsv'
    in_v2g = 'data/v2g/190201/evidence.tsv'
    outf = 'data/v2g/190201/v2g_evidence_gold_standards.merged.tsv'
    remove_ambiguous = True

    # Load gold standards
    gs = (
        pd.read_csv(in_gold_standards, sep='\t', header=0)
          .rename(columns={'gene_id': 'gold_standard_gene_id'})
    )
    print('Gold standard contains N variants: ', gs.shape[0])
    if remove_ambiguous:
        gs = gs.loc[~gs.varid_ambiguous, :]
        print('Removing ambiguous left N variants: ', gs.shape[0])

    # Load v2g
    v2g = pd.read_csv(in_v2g, sep='\t', header=0)

    # Calc normalised score
    v2g['normalisedOverallScore'] = calc_normalised_score(v2g)

    # Have to add missing gold standard genes to the v2g table
    gs_var_genes = ( gs.loc[:, ['varid', 'gold_standard_gene_id']]
                       .drop_duplicates() )
    v2g = pd.merge(v2g, gs_var_genes,
                   left_on=['varid', 'gene_id'],
                   right_on=['varid', 'gold_standard_gene_id'],
                   how='outer')
    v2g.loc[pd.isnull(v2g.gene_id), 'gene_id'] = (
        v2g.loc[pd.isnull(v2g.gene_id), 'gold_standard_gene_id'] )
    v2g = ( v2g.fillna(0)
               .drop('gold_standard_gene_id', axis=1) )

    # Merge gold standard and v2g
    data = pd.merge(gs, v2g, on='varid', how='left')
    data.to_csv(outf, sep='\t', index=None)

def calc_normalised_score(v2g):
    ''' Normalises the overall V2G score per varint by dividing by max per
        variant
    '''
    # Get max per variant
    max_p_var = (
        v2g.groupby('varid')
           .overallScore
           .max()
           .reset_index()
           .rename(columns={'overallScore': 'maxOverallScore'})
    )
    v2g = pd.merge(v2g, max_p_var, how='left')

    return v2g.overallScore / v2g.maxOverallScore

if __name__ == '__main__':

    main()
