#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import pandas as pd
import numpy as np
import seaborn as sns; sns.set()
import matplotlib
# matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def main():

    # Args
    in_data = 'data/v2g/190201/v2g_evidence_gold_standards.merged.tsv'
    out_dir = 'results/v2g/190201'

    os.makedirs(out_dir, exist_ok=True)

    # Load data
    data = pd.read_csv(in_data, sep='\t', header=0)

    # Only keep rows for gold standard evidence
    data = data.loc[data.gold_standard_gene_id == data.gene_id, :]

    # Make plot
    out_plot = os.path.join(out_dir, 'cumulative.plot.png')
    plot_cumulative_curve(data.normalisedOverallScore, out_plot)

    return 0

def plot_cumulative_curve(s, outf):
    ''' Plot score on x-axis and % of gold standards achieving score on Y-axis
    Params:
        s (pd.series): series of scores
        outf (file): file to output plot to
    '''
    # Create data
    x_list = np.linspace(0, 1, 21)
    y_list = [ ((s >= x).sum() * 100 / s.size) for x in x_list ]

    # Make plot
    ax = sns.lineplot(x=x_list, y=y_list)

    # Increase tick frequency
    ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())

    # Style
    sns.set(style='whitegrid')
    ax.set(
        xlabel='Normalised score (score / max(score))',
        ylabel='% Gold-standards ≥ Normalised score',
        title='V2G gold-standard validation'
    )

    # Save
    ax.get_figure().savefig(outf, dpi=150)

    return 0

if __name__ == '__main__':

    main()
