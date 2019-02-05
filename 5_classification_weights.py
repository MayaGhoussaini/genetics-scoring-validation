#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Ed Mountjoy
#

import sys
import os
import pandas as pd
import numpy as np
from functools import partial

# Classifiers
from sklearn.model_selection import GroupKFold
from sklearn.svm import LinearSVC
from sklearn.svm import SVC
from sklearn.linear_model import SGDClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neural_network import MLPClassifier
# Scorers
from sklearn.metrics import make_scorer, recall_score
from sklearn.model_selection import cross_val_score
# Stats
from statsmodels.discrete.discrete_model import Logit
from statsmodels.tools.tools import add_constant
# Scaler
from scripts.ml import maxabs_scaler

import sklearn
print('sklearn version: {}\n'.format(sklearn.__version__))

def main():

    # Args
    in_data = 'data/v2g/190201/v2g_evidence_gold_standards.merged.tsv'
    out_dir = 'results/v2g/190201'
    features = ['dhscor_thurman2012',
                'eqtl_gtex_v7',
                'fantom5_andersson2014',
                'fpred_vep',
                'pchic_javierre2016',
                'pqtl_sun2018']

    os.makedirs(out_dir, exist_ok=True)

    # Load data
    data = pd.read_csv(in_data, sep='\t', header=0)

    # Extract features and outcomes
    grps, X, y = extract_features(data, features, 'varid',
                                  remove_low_variance=True,
                                  group_scale=True)

    # Test classifiers
    # test_classifiers(X, y, grps, metric='balanced_accuracy', folds=5)
    # test_classifiers(X, y, grps, metric='roc_auc', folds=5)
    # test_classifiers(X, y, grps, metric='f1', folds=5)
    # test_classifiers(X, y, grps, metric='recall', folds=5)
    # recall_weighted = make_scorer(recall_score, average='weighted')
    # test_classifiers(X, y, grps, metric=recall_weighted, folds=5)

    # Learn feature weights
    learn_logistic_regression_coef(X, y, features)

    # Write example dataset
    # data['y'] = (data.gold_standard_gene_id == data.gene_id)
    # data.to_csv('data.test.tsv', sep='\t', index=None)

    return 0

def learn_logistic_regression_coef(X, y, x_labels):
    ''' Fit logistic regression and return feature coeffecieints
    '''
    X_const = pd.DataFrame(add_constant(X), columns=['Intercept'] + x_labels)
    clf = Logit(y, X_const)
    # Non regularised
    res = clf.fit()
    print('\n', res.summary())
    # # L1 regularised
    # res = clf.fit_regularized(method='l1')
    # print('\n', res.summary())
    # # l1_cvxopt_cp regularised
    # res = clf.fit_regularized(method='l1_cvxopt_cp')
    # print('\n', res.summary())

    return 0

def test_classifiers(X, y, grps, metric='balanced_accuracy', folds=5):
    ''' Use sklearn to test different classifiers

        metric choices: https://scikit-learn.org/stable/modules/model_evaluation.html#common-cases-predefined-values

    '''

    print('\nmetric=={0}'.format(metric))

    classifier_dict = {
        'LinearSVC': LinearSVC(),
        'SVC': SVC(gamma=2, C=1),
        'DecisionTreeClassifier': DecisionTreeClassifier(max_depth=5),
        'LogisticRegression_L2': LogisticRegression(penalty='l2', solver='lbfgs'),
        'LogisticRegression_L1': LogisticRegression(penalty='l1', solver='saga'),
        'KNeighborsClassifier': KNeighborsClassifier(),
        'RandomForestClassifier': RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
        'MLPClassifier': MLPClassifier(alpha=1),
        'AdaBoostClassifier': AdaBoostClassifier(),
        'SGDClassifier': SGDClassifier(max_iter=100, tol=1e-3)
    }

    # Create scores via cross-validation
    group_kfold = GroupKFold(n_splits=folds)
    for name, clf in classifier_dict.items():
        res = cross_val_score(clf, X, y,
                              scoring=metric,
                              groups=grps,
                              cv=group_kfold)
        print('{0}: {1:.3f} ± {2:.3f}'.format(name, np.mean(res), np.std(res)))

    return 0

def extract_features(df, feature_cols, group_col, remove_low_variance=True,
                     group_scale=False):
    ''' Returns feature (X), outcomes (y) and varid groups
    Params:
        df (pd.df)
        feature_cols (list)
        group_col (str)
        remove_low_variance (bool): remove rows with no variance in the predictors
        group_scale (bool):
    Returns
        groups, X, y (np.array)
    '''
    y = ( (df.gene_id == df.gold_standard_gene_id)
           .replace({True: 1, False: 0}) )
    grps = df[group_col]
    # Group normalise
    if not group_scale:
        X = df.loc[:, feature_cols]
    else:
        # Scale
        X = (
            df.loc[:, [group_col] + feature_cols]
              .groupby(group_col)
              .transform(maxabs_scaler)
        )
        # Fill NAs
        X = X.fillna(0)

    # Remove rows with no variance
    if remove_low_variance:
        to_keep = (X.var(axis=1) != 0)
        print('Removed {} samples (out of {}) with no feature variance'.format(
            (~to_keep).sum(), y.size))
        X = X.loc[to_keep, :]
        y = y.loc[to_keep]
        grps = grps.loc[to_keep]

    return grps.to_numpy(), X.to_numpy(), y.to_numpy()

if __name__ == '__main__':

    main()
