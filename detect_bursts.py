# -*- coding: utf-8 -*-

from __future__ import division

import argparse

import os
import sys

import numpy as np
import scipy.stats

import matplotlib.pyplot as plt

from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import FloatVector

try:
    import seaborn as sns
except:
    print 'using oldschool matplotlib'

import neo.io

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


importr('burstHSMM')
ro.r(get_probs_func_str)


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
    return 1.*len(spikes)/(spikes[-1] - spikes[0])


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


def detect_with_hist(spikes, args):
    spikes = np.array(spikes)

    t = 1./calc_discharge_rate(spikes)

    bins = np.arange(spikes[0],  spikes[-1]+t, t)
    hist_counts = np.bincount(np.histogram(spikes, bins)[0])
    count_values = list(range(len(hist_counts)))

    ch2_test = scipy.stats.chisquare(hist_counts, SPIKES_RV.pmf(count_values)*len(spikes))[1]
    skew_test = scipy.stats.skew(hist_counts)

    hist_counts = median_of_three_smoothing(hist_counts)

    if(ch2_test < 0.05 and skew_test > args.skewness):
        d = find_threshold_density(hist_counts)

        isi_t = t/d
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


def main(args):
    data_file = args.data_file

    if args.algorithm == 'hist':
        detect_func = detect_with_hist
    elif args.algorithm == 'hsmm':
        detect_func = detect_with_hsmm
    else:
        raise 'Unkown detect algorithm!'

    r = neo.io.NeuroExplorerIO(filename=data_file)
    blks = r.read(cascade=True, lazy=False)
    for blk in blks:
        for seg in blk.segments:
            for st in seg.spiketrains:
                spikes = np.array(st)
                if len(spikes) > 50:
                    burst_idx = detect_func(spikes, args)

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
    parser.add_argument('--min_spike_count', type=int, default=3,
                        help='Min spike count for burst bunch')

    args = parser.parse_args()

    main(args)
