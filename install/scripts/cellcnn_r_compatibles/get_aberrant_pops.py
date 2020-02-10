
"""
This code was modified by David Novak in 2019 for compatibility with tidycell via reticulate!

Get the original package here: https://imsb.ethz.ch/research/claassen/Software/cellcnn.html

Copyright 2016-2017 ETH Zurich, Eirini Arvaniti and Manfred Claassen.

"""

import copy
import os
import sys
import cPickle as pickle
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold, KFold
from cellCnn.utils import get_selected_cells
from cellCnn.plotting import plot_results, plot_filters, discriminative_filters
from cellCnn.model import CellCnn
from collections import Counter
from scipy import stats
from sklearn.cluster import DBSCAN
from sklearn.manifold import TSNE
from sklearn.preprocessing import MinMaxScaler
from sklearn.neighbors import NearestNeighbors
import seaborn as sns
import statsmodels.api as sm

def get_aberrant_pops(results_pointer,
                      fcs,
                      samples_pointer,
                      marker_names,
                      filter_diff_thres,
                      filter_response_thres,
                      select_discriminative_filters,
                      stat_test,
                      group_a,
                      group_b,
                      group_names,
                      regression,
                      train_samples,
                      train_phenotypes,
                      valid_samples,
                      valid_phenotypes):
    results = copy.deepcopy(results_pointer)
    samples = copy.deepcopy(samples_pointer)
    filter_info = get_filter_info(results, train_samples, train_phenotypes,
                                  marker_names,
                                  filter_diff_thres = filter_diff_thres,
                                  filter_response_thres = filter_response_thres,
                                  select_discriminative_filters = select_discriminative_filters,
                                  stat_test = stat_test,
                                  group_a = group_a, group_b = group_b,
                                  group_names = group_names,
                                  regression = regression)   

    nfilter = len(filter_info)
    fcs_info = np.array(pd.read_csv(fcs, sep=','))
    sample_names = [name.split('.fcs')[0] for name in list(fcs_info[:, 0])]
    # for each sample
    out = []
    for x, x_name in zip(samples, sample_names):
        flags = np.zeros((x.shape[0], 2*nfilter))
        columns = []
        # for each filter
        for i, (filter_idx, thres) in enumerate(filter_info):
            flags[:, 2*i:2*(i+1)] = get_selected_cells(
                results['selected_filters'][filter_idx], x, results['scaler'], thres, True)
            columns += ['filter_%d_continuous' % filter_idx, 'filter_%d_binary' % filter_idx]
        df = pd.DataFrame(flags, columns = columns)
        out.append(df)
    return out

def get_filter_info(results, samples, phenotypes, labels,
                filter_diff_thres, filter_response_thres,
                select_discriminative_filters,
                response_grad_cutoff=None,
                stat_test=None, log_yscale=False,
                group_a='group A', group_b='group B', group_names=None,
                regression=False, add_filter_response=False,
                percentage_drop_cluster=.1, min_cluster_freq=0.2):

    # get discriminative filter indices in consensus matrix
    print 'Fetching filter responses...'
    filters = results['selected_filters']

    # select the discriminative filters based on the validation set
    filter_diff = results['filter_diff']
    filter_diff[np.isnan(filter_diff)] = -1

    sorted_idx = []
    if select_discriminative_filters:
        sorted_idx = np.argsort(filter_diff)[::-1]
        filter_diff = filter_diff[sorted_idx]
        keep_idx = [sorted_idx[0]]
        for i in range(0, len(filter_diff)-1):
            if (filter_diff[i] - filter_diff[i+1]) < filter_diff_thres * filter_diff[i]:
                keep_idx.append(sorted_idx[i+1])
            else:
                break
    else:
        keep_idx = range(0, filters.shape[0])

    return_filters = []
    for i_filter in keep_idx:
        return_filters.append((i_filter, filter_response_thres))

    return return_filters
