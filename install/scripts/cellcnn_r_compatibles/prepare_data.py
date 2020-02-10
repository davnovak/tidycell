
"""
This code was modified by David Novak in 2019 for compatibility with tidycell via reticulate!

Get the original package here: https://imsb.ethz.ch/research/claassen/Software/cellcnn.html

Copyright 2016-2017 ETH Zurich, Eirini Arvaniti and Manfred Claassen.

"""

import os
import sys
import cPickle as pickle
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold, KFold
from cellCnn.utils import get_data, save_results, mkdir_p, get_selected_cells
from cellCnn.plotting import plot_results, plot_filters, discriminative_filters
from cellCnn.model import CellCnn

def prepare_data(fcs,
                 marker_names,
                 quant_normed,
                 arcsinh,
                 cofactor,
                 indir,
                 seed,
                 train_perc,
                 regression,
                 per_sample,
                 n_splits = None):
    # read in the data
    fcs_info = np.array(pd.read_csv(fcs, sep=','))

    # if the samples have already been pre-processed via quantile normalization
    # we should not perform arcsinh transformation
    if quant_normed:
        arcsinh = False

    samples, phenotypes = get_data(indir, fcs_info, marker_names,
                                    arcsinh, cofactor)

    # generate training/validation sets
    np.random.seed(seed)
    val_perc = 1 - train_perc
    if n_splits is None:
        n_splits = int(1. / val_perc)
    # stratified CV for classification problems
    if not regression:
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True)
    # simple CV for regression problems
    else:
        skf = KFold(n_splits=n_splits, shuffle=True)

    train, val = next(skf.split(np.zeros((len(phenotypes), 1)), phenotypes))
    #train = range(len(phenotypes))
    #val = None

    train_samples = [samples[i] for i in train]
    #valid_samples = [samples[i] for i in val] if val is not None else None
    valid_samples = [samples[i] for i in val]
    train_phenotypes = [phenotypes[i] for i in train]
    #valid_phenotypes = [phenotypes[i] for i in val] if val is not None else None
    valid_phenotypes = [phenotypes[i] for i in val]
# 
#     print '\nSamples used for model training:'
#     for i in train:
#         print fcs_info[i]
#     print '\nSamples used for validation:'
#     for i in val:
#         print fcs_info[i]
#     print

    # always generate multi-cell inputs on a per-sample basis for regression
    if regression:
        per_sample = True

    return train_samples, train_phenotypes, valid_samples, valid_phenotypes, marker_names, fcs_info, samples
