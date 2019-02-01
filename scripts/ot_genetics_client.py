#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

from gql import gql, Client
from gql.transport.requests import RequestsHTTPTransport

class OT_Genetics:
    ''' GraphQL client for Open Targets Genetics
    '''
    def __init__(self, api_endpoint=None, chunk_size=50):
        self.chunk_size = chunk_size
        if api_endpoint:
            url = api_endpoint
        # Setup graphql client
        _transport = RequestsHTTPTransport(
            url = 'https://genetics-api.opentargets.io/graphql',
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
