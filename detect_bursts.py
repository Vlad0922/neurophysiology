# -*- coding: utf-8 -*-

# HSMM based method described in paper "Detection of bursts in extracellular spike trains using hidden semi-Markov point process models"
# original source code for this method found here: http://www2.stat.duke.edu/~st118/Software/burstHSMM/

from __future__ import division, print_function

import argparse

import os

import numpy as np
import scipy.stats
import pandas as pd

from collections import Counter

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import FloatVector
from rpy2.robjects import pandas2ri
pandas2ri.activate()

try:
    import seaborn as sns
except:
    print('using oldschool matplotlib')

import neo.io

# importr('pracma')
# importr('sjemea')
# ro.r('source("3rdparty/CMA_method.R")')
# ro.r('source("3rdparty/PS_method.R")')

SPIKES_RV = scipy.stats.poisson(1.)

def plot_spikes(spikes, burst_spikes):
    f, axarr = plt.subplots(3, sharex=True)

    for i in range(spikes.shape[0]):
        axarr[0].plot([spikes[i], spikes[i]], [-1, 1], color='black', lw=0.5)

        if burst_spikes[i] == True:
            axarr[1].plot([spikes[i], spikes[i]], [-1, 1], color='red', lw=0.5)
            axarr[2].plot([spikes[i], spikes[i]], [-1, 1], color='red', lw=0.5)
        else:
            axarr[2].plot([spikes[i], spikes[i]], [-1, 1], color='black', lw=0.5)

    axarr[0].set_ylim([-5, 5])
    axarr[1].set_ylim([-5, 5])
    axarr[2].set_ylim([-5, 5])

    plt.show()

def clever_split(arr, step):
    res = list()
    
    idx = 0
    for l in np.arange(0, arr[~0], step):
        curr = list()
        while(idx < len(arr) and arr[idx] < l + step):
            curr.append(arr[idx])
            idx += 1
        
        res.append(curr)
    
    return res


def spiketrains_iterator(handler):
    for blk in handler.read(cascade=True, lazy=False):
        for seg in blk.segments:
            for st in seg.spiketrains:
                yield st.name, np.array(st)
                
def events_iterator(handler):
    for blk in handler.read(cascade=True, lazy=False):
        for seg in blk.segments:
            for st in seg.events:
                yield st.annotations['channel_name'], np.array(st)


def calc_discharge_rate(spikes):
    return 1.*len(spikes)/(spikes[~0] - spikes[0])

def bin_by_discharge(spikes):
    return 1./calc_discharge_rate(spikes)


def median_of_three_smoothing(hist):
    res = np.array(hist)
    res[0] = np.median([hist[0], hist[0], hist[1]])
    for i in range(1, hist.shape[0]-1):
        res[i] = np.median(hist[i-1:i+2])
    res[-1] = np.median([hist[-2], hist[-1], hist[-1]])

    return res


def find_by_bin(hist):
    d = None
    for i in range(1, hist.shape[0] - 1):
        if hist[i] < hist[i-1] and hist[i] <= hist[i+1]:
            d = i
            break

    return d


def find_by_slope(hist):
    d = None

    prev_slope = hist[1] - hist[0]
    for i in range(1, hist.shape[0]):
        curr_slope = hist[i] - hist[i-1]

        if abs(curr_slope) < (prev_slope) or (prev_slope < 0 and curr_slope > 0):
            d = i
            break

        prev_slope = curr_slope


    return d


def find_threshold_density(hist):
    d = find_by_bin(hist)

    if d is None:
        d = find_by_slope(hist)

    if d is None:
        raise 'Wow, cant find density threshold!'

    return d


def calc_bin_isi_std(spikes):
    isi = spikes[1:] - spikes[:-1]
    return np.std(isi)


def calc_bin_median(spikes):
    isi = spikes[1:] - spikes[:-1]
    return np.median(isi)


def detect_with_vitek(spikes, args, with_hist=False):
    try:
        coeff = 1.
        min_spikes = 3
        bin_func = bin_by_discharge

        spikes = np.array(spikes)
            
        t = bin_func(spikes)*coeff

        counts = list(map(len, clever_split(spikes, t)))
        nums, vals = zip(*sorted(Counter(counts).items()))

        vals = median_of_three_smoothing(np.array(vals, dtype=float))

        ch2_test = scipy.stats.chisquare(vals, SPIKES_RV.pmf(nums)*len(spikes))[1]
        skew_test = scipy.stats.skew(counts)

        burst_spikes = np.zeros(spikes.shape[0], dtype=bool)
        burst_bunches = list()
        burst_lens = list()  
        
        if(ch2_test < 0.05 and skew_test > 0.5):
            d = find_threshold_density(vals)

            isi_t = t/d
            
            isi = np.ediff1d(spikes)
            burst_isi = np.array(isi <= isi_t, dtype=bool)    
            
            prev = 0
            counter = 0
            for idx, i in enumerate(burst_isi):
                if i:
                    counter += 1
                else:
                    if counter >= min_spikes - 1:
                        burst_spikes[prev:idx+1] = True  
                        burst_bunches.append(spikes[prev:idx+1])
                        burst_lens.append(spikes[idx] - spikes[prev])

                    counter = 0   
                    prev = idx + 1
            
            if counter >= min_spikes - 1:  
                burst_spikes[prev:idx+1] = True  
                burst_bunches.append(spikes[prev:idx+1])
                burst_lens.append(spikes[idx] - spikes[prev])
            
        if with_hist:
            return burst_spikes, burst_bunches, burst_lens, vals
        else:
            return burst_spikes, burst_bunches, burst_lens
    except:
        return np.zeros(spikes.shape[0], dtype=int), [], []


def detect_plot_vitek(spikes, fname, args):
    bin_func_name = args.bin_func

    if bin_func_name == 'discharge':
        bin_func = bin_by_discharge
    elif bin_func_name == 'std':
        bin_func = calc_bin_isi_std
    elif bin_func_name == 'median':
        bin_func = calc_bin_median
    else:
        raise 'Unknown bin function!'

    st = np.array(spikes)
    bin_size = bin_func(st)*args.bin_coeff

    burst_spikes, burst_bunches, _, hist = detect_with_vitek(st, args=args, with_hist=True)

    d = find_threshold_density(hist)
    isi_t = bin_size/d

    print(d)
    print(hist)

    isi = st[1:] - st[:-1]
    
    fig = plt.figure(figsize=(15,10))
    gridspec.GridSpec(3,3)

    plt.subplot2grid((3, 3), (0,1))
    plt.bar(np.arange(hist.shape[0]), hist)
    plt.plot(SPIKES_RV.pmf(np.arange(hist.shape[0]))*len(st), color='black')
    plt.title('SDH, method: {}, bin size: {}'.format(args.bin_func, round(bin_size, 3)))
    
    plt.subplot2grid((3, 3), (1,0), colspan=3)
    for s in st:
        plt.plot([s, s], [-1, 1], color='black', lw=0.3)
        
    plt.ylim([-5, 5])      
    plt.xlim([-0.5, st[~0] + 0.5])
    plt.title('raw spikes, name: {}'.format(fname))
    
    plt.subplot2grid((3, 3), (2,0), colspan=3)
    for s in st[burst_spikes]:
        plt.plot([s, s], [-1, 1], color='red', lw=0.3)
        
    for b in burst_bunches:
        plt.plot([b[0], b[~0]], [-1.2, -1.2], color='green', lw=0.5)
    
    plt.ylim([-5, 5])
    plt.xlim([-0.5, st[~0] + 0.5])
    plt.title('burst spikes, ai: {}, burst isi: {}'.format(round(np.median(isi)/np.mean(isi), 3), round(isi_t, 3)))  
    plt.show()  



def detect_with_hsmm(spikes, args):
    spikes = np.array(spikes)
    isi = np.ediff1d(spikes)
    isi_r = ro.FloatVector(isi)

    ro.r.assign('isi', isi_r)
    ro.r('thrash <- capture.output(model <- fit.hsmm(isi))')

    probs = np.array(ro.r('get_probs_hsmm(model)'))

    return find_burst_spikes(spikes, np.array(probs > args.prob_threshold), args.min_spike_count)


def detect_with_rs(spikes, args):
    spikes = np.array(spikes)
    start, length, RS = rank_surprise.burst(spikes, limit=50e-3, RSalpha=1.)

    burst_spikes = np.zeros(spikes.shape[0])
    for s, l in zip(start, length):
        e = s + l
        burst_spikes[s:e] = True

    return burst_spikes


def detect_with_logisi(spikes, args):
    spikes_r = ro.FloatVector(spikes)

    ro.r.assign('spikes', spikes_r)
    res = pd.DataFrame(ro.r('data.frame(logisi.pasq.method(spikes, cutoff={}))'.format(args.logisi_cutoff)))

    burst_spikes = np.zeros(spikes.shape[0])

    try:
        idx = zip(map(int, res['beg']), map(int, res['end']))

        for s, e in idx:
            s -= 1
            burst_spikes[s:e] = True

        burst_bunches = list()
        for s, e in idx:
            s -= 1
            burst_bunches.append(np.array(spikes[s:e]))

        return burst_spikes, burst_bunches
    except:
        return burst_spikes, []


def detect_with_ps(spikes, args):
    try:
        thresh = args.si_thresh
    except:
        thresh = 3

    spikes_r = ro.FloatVector(spikes)

    ro.r.assign('spikes', spikes_r)
    res = pd.DataFrame(ro.r('data.frame(PS.method(spikes, si.thresh={}))'.format(thresh)))

    burst_spikes = np.zeros(spikes.shape[0])

    try:
        idx = list(zip(map(int, res['beg']), map(int, res['end'])))
        burst_lens = list()

        for s, e in idx:
            s -= 1
            burst_spikes[s:e] = True
            burst_lens.append(spikes[e-1] - spikes[s])

        burst_bunches = list()
        for s, e in idx:
            s -= 1
            burst_bunches.append(np.array(spikes[s:e]))

        return burst_spikes, burst_bunches, burst_lens
    except:
        return burst_spikes, [], []


def detect_with_cma(spikes, args):
    spikes_r = ro.FloatVector(spikes)

    ro.r.assign('spikes', spikes_r)
    res = pd.DataFrame(ro.r('data.frame(CMA.method(spikes))'))

    burst_spikes = np.zeros(spikes.shape[0], dtype=int )

    try:
        idx = list(zip(map(int, res['beg']), map(int, res['end'])))
        burst_lens = list()
        burst_bunches = list()

        for s, e in idx:
            s -= 1

            burst_spikes[s:e] = 1
            burst_lens.append(spikes[e-1] - spikes[s])
            burst_bunches.append(np.array(spikes[s:e]))

        return burst_spikes, burst_bunches, burst_lens
    except:
        return burst_spikes, [], []


def main(args):
    data_file = args.data_file
    algo = args.algorithm

    if algo == 'vitek':
        detect_func = detect_plot_vitek
    elif algo == 'hsmm':
        detect_func = detect_with_hsmm
    elif algo == 'rs':
        detect_func = detect_with_rs
    elif algo == 'logisi':
        detect_func = detect_with_logisi
    else:
        raise 'Unkown detect algorithm!'

    ext = data_file[-3:].lower()
    if ext == 'nex':
        r = neo.io.NeuroExplorerIO(filename=data_file)
    elif ext == 'smr':
        r = neo.io.Spike2IO(filename=data_file)

    for name, st in events_iterator(r):
        spikes = np.array(st)
        if len(spikes) > 50:
            detect_func(spikes, name, args)


    for name, st in spiketrains_iterator(r):
        spikes = np.array(st)
        if len(spikes) > 50:
            detect_func(spikes, name, args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Detect bursts in spike data')

    parser.add_argument('--data_file', type=str, required=True,
                        help='Nex data file with spikes')
    parser.add_argument('--algorithm', type=str, default='vitek',
                        help='Which algorithm use in spike detection')
    parser.add_argument('--skewness', type=float, default=1.,
                        help='Skewness param')
    parser.add_argument('--prob_threshold', type=float, default=0.5,
                        help='probability threshold for hsmm algorithm')
    parser.add_argument('--plot', action='store_true', default=True,
                        help='Plot burst and non burst spikes?')
    parser.add_argument('--min_spikes', type=int, default=3,
                        help='Min spike count for burst bunch')
    parser.add_argument('--bin_func', type=str, default='discharge',
                        help='Function for bin calculation')
    parser.add_argument('--bin_coeff', type=float, default=1.,
                        help='Coefficient of disharge rate multiplication')

    args = parser.parse_args()

    main(args)
