#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

def maxabs_scaler(x):
    ''' Max absolute scaler
    Params:
        x (pd.Series)
    Returns
        pd.Series
    '''
    return (x / float(x.abs().max()))
