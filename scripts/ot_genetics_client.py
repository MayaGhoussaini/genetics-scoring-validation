#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

from gql import gql, Client
from gql.transport.requests import RequestsHTTPTransport
import pandas as pd

class OT_Genetics:
    ''' GraphQL client for Open Targets Genetics
    '''
    def __init__(self, api_endpoint=None, chunk_size=50):
        self.chunk_size = chunk_size
        if api_endpoint:
            otg_url = api_endpoint
        else:
            otg_url = 'https://genetics-api.opentargets.io/graphql'

        # Setup graphql client
        _transport = RequestsHTTPTransport(
            url = otg_url,
            use_json = True
        )
        
        self.client = Client(
            transport=_transport,
            fetch_schema_from_transport=True,
        )

    def rsid_to_varid(self, rsid_list):
        ''' Convert rsids to a dictionary {rsid -> [varid1, varid2, ...]}
        '''
        # Make query
        res = {}
        for chunk in get_chunk(iter(rsid_list), self.chunk_size):
            # Build inner part of query
            inner = '\n'.join([
                '''{rsid}:
                search(queryString:"{rsid}") {{
                     variants {{
                        variant {{
                            id
                        }}
                    }}
                }}
                    '''.format(rsid=rsid)
                for rsid in chunk ])
            # Build query
            query = 'query search_rsid{ ' + inner + '}'
            # Query and return
            chunk_res = self.client.execute(gql(query))
            res.update(chunk_res)

        # Convert to dictionary
        rsid_map = {}
        for rsid in res:
            rsid_map[rsid] = [
                record['variant']['id'] for record in res[rsid]['variants']
            ]

        return rsid_map

    def get_v2g_evidence_for_varids(self, varid_list):
        ''' Gets v2g evidence for a list of variant IDs
        Returns:
            pd.df
        '''
        # Make query
        res = {}
        for chunk in get_chunk(iter(varid_list), self.chunk_size):
            # Build inner part of query
            inner = '\n'.join([
                '''v{varid}:
                genesForVariant(variantId:"{varid}") {{
                     gene {{
                       id
                     }}
                     overallScore
                     qtls {{
                       typeId
                       sourceId
                       aggregatedScore
                     }}
                     intervals {{
                       typeId
                       sourceId
                       aggregatedScore
                     }}
                     functionalPredictions {{
                       typeId
                       sourceId
                       aggregatedScore
                     }}
                }}
                    '''.format(varid=varid)
                for varid in chunk ])
            # Build query
            query = 'query get_v2g{ ' + inner + '}'
            # Query and return
            chunk_res = self.client.execute(gql(query))
            res.update(chunk_res)

        return v2g_to_df(res)

def v2g_to_df(res):
    ''' Convert the response from V2G query to a pandas df
    '''
    rows = []
    for varid, records in res.items():
        varid = varid.lstrip('v')
        for record in records:
            row = {'varid': varid}
            row['gene_id'] = record['gene']['id']
            row['overallScore'] = record['overallScore']
            for type in ['qtls', 'intervals', 'functionalPredictions']:
                for entry in record[type]:
                    key = '_'.join([entry['typeId'], entry['sourceId']])
                    row[key] = entry['aggregatedScore']
            rows.append(row)
    df = pd.DataFrame.from_dict(rows)

    # Move columns to front
    df = move_cols_to_front(df, ['varid', 'gene_id', 'overallScore'])

    return df

def move_cols_to_front(df, cols):
    ''' Moves a list of columns to the front of a pandas df
    '''
    existing_columns = df.columns.tolist()
    for col in cols:
        existing_columns.remove(col)

    return df.loc[:, cols + existing_columns]

def get_chunk(it, n):
    ''' Returns chunks of an iterator
    '''
    try:
        while True:
            xs = []  # The buffer to hold the next n items
            for _ in range(n):
                xs.append(next(it))
            yield xs
    except StopIteration:
        yield xs
