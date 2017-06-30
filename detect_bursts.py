# -*- coding: utf-8 -*-

from __future__ import division

import os
import sys

import numpy as np
import scipy.stats

import matplotlib.pyplot as plt

try:
    import seaborn as sns
except:
    print 'using oldschool matplotlib'

import neo.io


SPIKES_RV = scipy.stats.poisson(1.)


def calc_discharge_rate(spikes):
    return 1.*len(spikes)/(spikes[-1] - spikes[0])

def median_of_three_smoothing(hist):
    res = np.array(hist)
    res[0] = np.median([hist[0], hist[0], hist[1]])
    for i in range(1, hist.shape[0]-1):
        res[i] = np.median(hist[i-1:i+2])
    res[-1] = np.median([hist[-2], hist[-1], hist[-1]])

    return res


def find_burst_spikes(spikes, isi, isi_t):
    burst_spikes = np.zeros(spikes.shape[0], dtype=bool)

    prev = 0
    counter = 0
    for i in range(isi.shape[0]):
        if isi[i] <= isi_t:
            counter += 1
        else:
            if counter > 2:
                burst_spikes[prev:i+1] = True
            counter = 0
            prev = i+1

    if counter > 2:
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


def plot_discharge_hist(spikes, title=None, skewness=1., plot=False):
    spikes = np.array(spikes)

    t = 1./calc_discharge_rate(spikes)

    bins = np.arange(spikes[0],  spikes[-1]+t, t)
    hist_counts = np.bincount(np.histogram(spikes, bins)[0])
    count_values = list(range(len(hist_counts)))

    ch2_test = scipy.stats.chisquare(hist_counts, SPIKES_RV.pmf(count_values)*len(spikes))[1]
    skew_test = scipy.stats.skew(hist_counts)

    print '{}, ch2: {}, skew: {}'.format(title, ch2_test, skew_test)

    hist_counts = median_of_three_smoothing(hist_counts)

    if(ch2_test < 0.05 and skew_test > skewness):
        print 'Wow, we have bursts!'
        d = find_threshold_density(hist_counts)

        isi_t = t/d
        isi = np.ediff1d(spikes)

        burst_spikes = find_burst_spikes(spikes, isi, isi_t)

        print 'burst neurons count: {}/{}'.format(np.count_nonzero(burst_spikes), len(spikes))

        if plot:
            f, axarr = plt.subplots(2, sharex=True)

            for i in range(spikes.shape[0]):
                axarr[0].plot([spikes[i], spikes[i]], [-1, 1], color='black', lw=0.5)

                if burst_spikes[i] == True:
                    axarr[1].plot([spikes[i], spikes[i]], [-1, 1], color='red', lw=0.5)

            for i in range(isi.shape[0]):
                if isi[i] <= isi_t:
                    axarr[0].plot([spikes[i+1], spikes[i]], [0, 0], color='green', lw=1)
                    axarr[1].plot([spikes[i+1], spikes[i]], [0, 0], color='green', lw=1)

            axarr[0].set_ylim([-5, 5])
            axarr[1].set_ylim([-5, 5])

            plt.title('{}, ch2: {}, skew: {}'.format(title, ch2_test, skew_test))
            plt.show()

    if plot:
        plt.bar(count_values, hist_counts, color='green')
        if title is not None:
            plt.title('{}, ch2: {}, skew: {}'.format(title, ch2_test, skew_test))

        plt.show()


def main():
    data_file = sys.argv[1]
    skewness = float(sys.argv[2])

    if len(sys.argv) == 4:
        plot = (sys.argv[3] == 'plot')
    else:
        plot = False

    r = neo.io.NeuroExplorerIO(filename=data_file)
    blks = r.read(cascade=True, lazy=False)
    for blk in blks:
        for seg in blk.segments:
            for st in seg.spiketrains:
                spikes = np.array(st)
                if len(spikes) > 50 and not st.name.startswith('fon'):
                    plot_discharge_hist(spikes, st.name, skewness, plot)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'wrong usage!'
        print 'usage: python detect_bursts.py <nex_data_file> <skewness> <plot[optional]>'
        exit(1)

    main()
