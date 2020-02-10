
"""
This code was modified by David Novak in 2019 for compatibility with tidycell via reticulate!

Get the original package here: https://imsb.ethz.ch/research/claassen/Software/cellcnn.html

Copyright 2016-2017 ETH Zurich, Eirini Arvaniti and Manfred Claassen.

"""

import os
import sys
from collections import Counter
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import DBSCAN
from sklearn.manifold import TSNE
from sklearn.preprocessing import MinMaxScaler
from sklearn.neighbors import NearestNeighbors
from sklearn.utils import shuffle
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import ImageGrid
import seaborn as sns
from cellCnn.utils import mkdir_p
import statsmodels.api as sm
try:
    from cellCnn.utils import create_graph
except ImportError:
    pass


def filter_info(results, samples, phenotypes, labels,
                filter_diff_thres=.2, filter_response_thres=0, response_grad_cutoff=None,
                stat_test=None, log_yscale=False,
                group_a='group A', group_b='group B', group_names=None, tsne_ncell=10000,
                regression=False, add_filter_response=False,
                percentage_drop_cluster=.1, min_cluster_freq=0.2):
    # number of measured markers
    nmark = samples[0].shape[1]

    if results['selected_filters'] is not None:
        print 'Loading the weights of consensus filters.'
        filters = results['selected_filters']
    else:
        sys.exit('Consensus filters were not found.')

    # get discriminative filter indices in consensus matrix
    keep_idx = discriminative_filters(results, filter_diff_thres)

    # encode the sample and sample-phenotype for each cell
    sample_sizes = []
    per_cell_ids = []
    for i, x in enumerate(samples):
        sample_sizes.append(x.shape[0])
        per_cell_ids.append(i * np.ones(x.shape[0]))
    
    return_filters = []
    for i_filter in keep_idx:
        t = filter_response_thres
        return_filters.append((i_filter, t))
                
    print 'Done.\n'
    return return_filters


def discriminative_filters(results, filter_diff_thres):
    filters = results['selected_filters']
    # select the discriminative filters based on the validation set
    if 'filter_diff' in results:
        filter_diff = results['filter_diff']
        filter_diff[np.isnan(filter_diff)] = -1

        sorted_idx = np.argsort(filter_diff)[::-1]
        filter_diff = filter_diff[sorted_idx]
        keep_idx = [sorted_idx[0]]
        for i in range(0, len(filter_diff)-1):
            if (filter_diff[i] - filter_diff[i+1]) < filter_diff_thres * filter_diff[i]:
                keep_idx.append(sorted_idx[i+1])
            else:
                break
    elif 'filter_tau' in results:
        filter_diff = results['filter_tau']
        filter_diff[np.isnan(filter_diff)] = -1

        sorted_idx = np.argsort(filter_diff)[::-1]
        filter_diff = filter_diff[sorted_idx]
        keep_idx = [sorted_idx[0]]
        for i in range(0, len(filter_diff)-1):
            if (filter_diff[i] - filter_diff[i+1]) < filter_diff_thres * filter_diff[i]:
                keep_idx.append(sorted_idx[i+1])
            else:
                break
    # if no validation samples were provided, keep all consensus filters
    else:
        filters = results['selected_filters']
        keep_idx = range(filters.shape[0])
    return keep_idx
