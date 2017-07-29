# -*- coding: utf-8 -*-

# HSMM based method described in paper "Detection of bursts in extracellular spike trains using hidden semi-Markov point process models"
# original source code for this method found here: http://www2.stat.duke.edu/~st118/Software/burstHSMM/

from __future__ import division

import argparse

import os

import numpy as np
import scipy.stats
import pandas as pd

import matplotlib.pyplot as plt

from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import FloatVector
from rpy2.robjects import pandas2ri
pandas2ri.activate()

try:
    import seaborn as sns
except:
    print 'using oldschool matplotlib'

import neo.io

import rank_surprise

get_probs_func_str = 'get_probs_hsmm <- function(out)                           \n\
                        {                                                       \n\
                          subset <- 1:out$mcmc[2]                               \n\
                                                                                \n\
                          p.states <- 1 - apply(out$states[, subset], 1, mean)  \n\
                          pars <- out$pars[, subset]                            \n\
                                                                                \n\
                          if(median(pars[2,] - pars[4,]) > 0){                  \n\
                            p.states <- 1 - p.states                            \n\
                          }                                                     \n\
                                                                                \n\
                          return(p.states)                                      \n\
                        }'

importr('pracma')
ro.r(get_probs_func_str)
ro.r('source("logisi.R")')

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


def calc_discharge_rate(spikes):
    return 1.*len(spikes)/(spikes[~0] - spikes[0])


def median_of_three_smoothing(hist):
    res = np.array(hist)
    res[0] = np.median([hist[0], hist[0], hist[1]])
    for i in range(1, hist.shape[0]-1):
        res[i] = np.median(hist[i-1:i+2])
    res[-1] = np.median([hist[-2], hist[-1], hist[-1]])

    return res


def find_burst_spikes(spikes, burst_isi, min_count):
    min_count -= 1

    burst_spikes = np.zeros(spikes.shape[0], dtype=bool)

    prev = 0
    counter = 0
    for i in range(burst_isi.shape[0]):
        if burst_isi[i]:
            counter += 1
        else:
            if counter > min_count:
                burst_spikes[prev:i+1] = True
            counter = 0
            prev = i+1

    if counter > min_count:
        burst_spikes[prev:i+1] = True

    return burst_spikes


def find_by_bin(hist):
    d = None
    for i in range(1, hist.shape[0] - 1):
        if hist[i] < hist[i-1] and hist[i] <= hist[i+1]:
            d = i
            break

    return d


def find_by_slope(hist):
    d = None
    slope, _, _, _, _ = scipy.stats.linregress(list(range(len(hist))), hist)
    slope = abs(slope)
    for i in range(1, hist.shape[0]):
        if abs(hist[i] - hist[i-1]) < slope:
            d = i
            break

    return d


def find_threshold_density(hist):
    d = find_by_bin(hist)

    if d is None:
        d = find_by_slope(hist)

    if d is None:
        raise 'Wow, cant find density threshold!'

    return d


def find_burst_bunches(spikes, burst_mask):
    res = list()
    curr = list()
    for i in range(burst_mask.shape[0]):
        if burst_mask[i]:
            curr.append(spikes[i])
        else:
            if len(curr) != 0:
                res.append(curr)
            curr = list()
    
    return res


def detect_with_vitek(spikes, args):
    spikes = np.array(spikes)

    t = 1./calc_discharge_rate(spikes)

    bins = np.arange(spikes[0],  spikes[~0] + t, t)
    hist_counts = np.bincount(np.histogram(spikes, bins)[0])
    count_values = list(range(len(hist_counts)))

    hist_counts = median_of_three_smoothing(hist_counts)

    ch2_test = scipy.stats.chisquare(hist_counts, SPIKES_RV.pmf(count_values)*len(spikes))[1]
    skew_test = scipy.stats.skew(hist_counts)

    print ch2_test, skew_test

    if(ch2_test < 0.05 and skew_test > args.skewness):
        d = find_threshold_density(hist_counts)

        isi_t = t/(d-1)
        isi = np.ediff1d(spikes)
        burst_isi = np.array(isi <= isi_t)

        return find_burst_spikes(spikes, burst_isi, args.min_spike_count)

    return np.zeros(spikes.shape[0])


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
    res = pd.DataFrame(ro.r('data.frame(logisi.pasq.method(spikes))'))

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


def main(args):
    data_file = args.data_file
    algo = args.algorithm

    if algo == 'vitek':
        detect_func = detect_with_vitek
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
    blks = r.read(cascade=True, lazy=False)
    for blk in blks:
        for seg in blk.segments:
            for st in seg.spiketrains:
                spikes = np.array(st)
                if len(spikes) > 50:
                    burst_idx, bunches = detect_func(spikes, args)

                    if args.plot:
                        plot_spikes(spikes, burst_idx)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Detect bursts in spike data')

    parser.add_argument('--data_file', type=str, required=True,
                        help='Nex data file with spikes')
    parser.add_argument('--algorithm', type=str, default='hist',
                        help='Which algorithm use in spike detection')
    parser.add_argument('--skewness', type=float, default=1.,
                        help='Skewness param')
    parser.add_argument('--prob_threshold', type=float, default=0.5,
                        help='probability threshold for hsmm algorithm')
    parser.add_argument('--plot', action='store_true', default=False,
                        help='Plot burst and non burst spikes?')
    parser.add_argument('--min_spike_count', type=int, default=2,
                        help='Min spike count for burst bunch')

    args = parser.parse_args()

    main(args)
