# -*- coding: utf-8 -*-

import sys
import os

import nexfile

import numpy as np
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, set_link_color_palette

import scipy as sp

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches   

from collections import Counter

import pandas as pd

from sklearn.preprocessing import normalize
from sklearn.cluster import KMeans

from collections import defaultdict

from scipy.stats import variation

import re
import pickle

import argparse
import scipy.spatial.distance as ssd

from spike_filter import apply_intervals, get_intervals, get_spiketrains

PATTERN_COLORS = {'tonic': 'green', 'irregular': 'yellow', 'burst': 'orange', 'pause': 'red'}

def JSD(P, Q):
    _P = P / np.linalg.norm(P, ord=1)
    _Q = Q / np.linalg.norm(Q, ord=1)
    _M = 0.5 * (_P + _Q)
    return 0.5 * (sp.stats.entropy(_P, _M) + sp.stats.entropy(_Q, _M))


def merge_st(st_list):
    offset = st_list[0][~0]

    for i in range(1, len(st_list)):
        st_list[i] = st_list[i][1:] - st_list[i][0] + offset
        offset = st_list[i][~0]

    return np.concatenate(st_list)


def calc_ai(spikes):
    intervals = spikes[1:] - spikes[:-1]
    return np.median(intervals)/np.mean(intervals)


def calc_cv(spikes):
    intervals = spikes[1:] - spikes[:-1]
    return variation(intervals)


def extract_depth(filename):
    vals = re.findall(r"\d+\.\d+", filename.replace('_', '.'))
    return float(vals[0])

def get_spikes(dist_dir):
    all_data = defaultdict(list)
    for root, subdirs, files in os.walk(dist_dir):
        for full_name, f_name in [(os.path.join(root, f_name), f_name) for f_name in files]:
            patient = full_name.split(os.sep)[~1]

            ext = full_name.split('.')[~0].lower()
            depth = extract_depth(f_name)

            if not('nex' in ext.lower()):
                continue

            r = nexfile.Reader()
            file_data = r.ReadNexFile(filePath=full_name)

            spiketrains = list(get_spiketrains(file_data))
            intervals = list(get_intervals(file_data))

            for spiketrain_name, interval_name, spikes in apply_intervals(spiketrains, intervals, fixed_interval_name=args.interval_name):
                if len(spikes) > 50 and (spikes[~0] - spikes[0] > 5.):
                    df = dict()
                    df['spikes'] = spikes
                    df['patient'] = patient
                    df['data_name'] = spiketrain_name
                    df['doc_name'] = f_name
                    df['interval_name'] = interval_name
                    df['depth'] = depth
                    df['firing rate'] = 1.*len(spikes)/(spikes[~0] - spikes[0])
                    df['CV'] = calc_cv(spikes)
                    df['AI'] = calc_ai(spikes)
                    
                    for k in df:
                        all_data[k].append(df[k])

    return all_data


def get_sdh(st, norm = True):
    isimean = np.mean(st[1:] - st[:-1])
    rec = (st[~0] - st[0])/len(st)
    
    binsize = (isimean + rec)/2
    
    cnt, _ = np.histogram(st, bins=np.arange(st[0], st[~0] + binsize, binsize))
    counts = Counter(cnt)
    
    res = np.zeros(max(counts.keys()) + 1)
    
    for idx, val in counts.items():
        res[idx] = val
    
    if(norm):
        res /= sum(res)
    
    return res 

def pad_to_size(arr, sz, val):
    res = np.full(sz, val, dtype=float)
    res[:len(arr)] = arr[:sz]
    
    return res


def find_cut(Z, num, l = 0.05, r = 10.):
    t = (l + r)/2
    curr_cut = fcluster(Z, t, criterion='distance')
    curr_num = len(np.unique(curr_cut))
    
    while curr_num != num:        
        if curr_num > num:
            l = t
            t = (t + r)/2
        elif curr_num < num:
            r = t
            t = (t + l)/2
            
        t = (l + r)/2
        curr_cut = fcluster(Z, t, criterion='distance')
        curr_num = len(np.unique(fcluster(Z, t, criterion='distance')))
    
    return t, curr_cut


def plot_tree(Z_curr, t, pattern_names):
    colors = [PATTERN_COLORS[pname] for pname in pattern_names]
    set_link_color_palette(colors)

    fig_tree, axes_dend = plt.subplots(ncols = 2, figsize=(20, 10))
    _ = dendrogram(Z_curr, color_threshold=t, above_threshold_color="grey", ax=axes_dend[0])

    optim = np.argmax(Z_curr[:,2][::-1][:9] / Z_curr[:,2][::-1][1:10]) + 1   
    
    axes_dend[1].plot(np.arange(1, 11), Z_curr[:,2][::-1][:10])
    axes_dend[1].plot([optim + 1], [Z_curr[:, 2][::-1][optim]], marker='o', color='red', markersize=10)
    axes_dend[1].set_title('Cost of merging')
    axes_dend[1].set_xlabel('Number of clusters')
    axes_dend[1].set_ylabel('Sum of squares')
    axes_dend[1].set_xticks(np.arange(1, 11))
    axes_dend[1].set_xticklabels(np.arange(1, 11))

    axes_dend[0].set_title('Dendogram')
    axes_dend[0].set_ylabel('Sum of squares')
    axes_dend[0].set_xticks([])
    
    
    patches = [mpatches.Patch(color=PATTERN_COLORS[pname], label=pname) for pname in pattern_names]
    axes_dend[0].legend(handles=patches)

    fig_tree.savefig('dendogram.png')


def get_structure(d):
    if (d >= 0 and d <= 6):
        return 'GPi'
    else:
        return 'GPe'
    

def get_pat_name(l, pattern_names, cmap):
    return pattern_names[cmap.index(l)]


def plot_cv_ai(df):
    fig, ax = plt.subplots(figsize=(10,10))
    def _plot_single(grouped, pname):
        pcol = PATTERN_COLORS[pname]
        ax.scatter(grouped['CV'], grouped['AI'], label=pname, color=pcol)
    
    for (pname, group) in df.groupby('Pattern'):
        _plot_single(group, pname)
    
    ax.set_xlabel('CV')
    ax.set_ylabel('AI')
    ax.legend()

    fig.savefig('cv_ai.png')


def main(args):
    data = get_spikes(args.data_dir)

    all_sdh = np.array([pad_to_size(get_sdh(st), 20, 0.0) for st in data['spikes']])

    M = sp.spatial.distance.pdist(all_sdh, metric=JSD)
    M_dists = ssd.squareform(M)
    
    n_clusts = args.clusters
    Z_curr = linkage(normalize(all_sdh), 'ward')
    t = Z_curr[-(n_clusts-1),2]
    #t, hcut_curr = find_cut(Z_curr, n_clusts)
    hcut_curr = fcluster(Z_curr, t=n_clusts, criterion='maxclust') - 1
    
    mdl = KMeans(n_clusters=n_clusts).fit(all_sdh)
    cut_curr = mdl.labels_        

    h_centroids = [all_sdh[hcut_curr == l].mean(axis=0) for l in np.unique(hcut_curr)]
    hcentroids_indexes = [i[0] for i in sorted(enumerate(h_centroids), key=lambda x:x[1][0]/x[1][1])]

    centroids = sorted(mdl.cluster_centers_, key = lambda s: s[0]/s[1])
    for i in range(len(centroids)):
        centroids[i][centroids[i] < 0] = 0

    centroids_indexes = [i[0] for i in sorted(enumerate(mdl.cluster_centers_), key=lambda x:x[1][0]/x[1][1])]

    if len(centroids) == 2:
        pattern_names = ['tonic', 'burst']
    elif len(centroids) == 3:
        pattern_names = ['tonic', 'burst', 'pause',]
    else:
        pattern_names = ['tonic', 'irregular', 'burst', 'pause']
    
    del data['spikes']
    df = pd.DataFrame(data)

    df['Structure'] = df['depth'].apply(get_structure)
    df['Pattern'] = hcut_curr
    df['Pattern'] = df['Pattern'].apply(get_pat_name, pattern_names=pattern_names, cmap=hcentroids_indexes)
    
    df = df[['patient', 'doc_name', 'data_name', 'interval_name', 'depth', 'Structure', 'Pattern', 'AI', 'CV', 'firing rate']]
    
    target_fname = args.dist_file + '.xlsx'
    df.to_excel(target_fname, index=False, header=True)
    
    if args.plot:
        plot_tree(Z_curr, t, pattern_names)
        plot_cv_ai(df)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Detect bursts in spike data')

    parser.add_argument('--data_dir',   type=str, required=True, help='Input directory')
    parser.add_argument('--dist_file',  type=str, required=True, help='File with results')
    parser.add_argument('--clusters',   type=int, default=4, help='Number of clusters')
    parser.add_argument('--interval_name', type=str, default=None)
    parser.add_argument('--plot', help='File with results', action='store_true')

    args = parser.parse_args()

    main(args)