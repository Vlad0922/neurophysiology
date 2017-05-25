# -*- coding: utf-8 -*-

from __future__ import division


import numpy as np
from scipy.special import psi
from scipy.stats import variation, kurtosis, skew
from scipy.signal import periodogram

FILTER_TIME = 1
FILTER_NUMBER = 0

DEFAULT_BBH = 0.016
DEFAULT_BBL = 0.01
DEFAULT_WIN_SIZE = 13
DEFAULT_REPEAT = 3
DEFAULT_FILT_LOWER = -1
DEFAULT_FILT_UPPER = -1
DEFAULT_FILT_TYPE = FILTER_TIME
DEFAULT_PAUSE_INDEX_BOUND = 50
DEFAULT_MODALIRITY_BOUND = 10
DEFAULT_PAUSE_RATIO_BOUND = 50
DEFAULT_BURST_BEHAVIOUR_BOUND = 100
DEFAULT_STEP = 1
DEFAULT_FREQUENCY = 1000.


def idx_of_nearest(arr, val):
    arr = np.array(arr)
    val, idx = min((val, idx) for (idx, val) in enumerate(np.abs(arr-val)))
    return idx


def sec_to_timestamps(spikes, freq):
    return np.array((spikes-np.floor(spikes[0]))*freq, dtype=np.int32)


def sec_to_ms(sec):
    return sec*1000


def ms_to_sec(ms):
    return 1.*ms/1000


def filter_spikes(spike_data, begin, end):
    return spike_data[(spike_data >= begin) & (spike_data <= end)]
    

def calc_intervals(spike_data, step=1):
    return spike_data.flat[step:] - spike_data.flat[:-step]


def calc_cv(intervals):
    return variation(intervals)


def calc_nu(intervals, wsize=DEFAULT_WIN_SIZE):
    n = len(intervals)
    m = wsize

    if m > n:
        m = n//2

    int_mean = np.mean(intervals)
    int_sorted = np.sort(intervals)

    bi2_vas_int = np.fromfunction(lambda x: psi(x+m), (m,))
    bi2_vas = np.log(2*m/n) - (1-2*m/n)*psi(2*m)+psi(n+1)-(2/n)*sum(bi2_vas_int)

    tadd = np.zeros((n,))
    tadd[n-m:] = int_sorted[-1]
    tadd[:n-m] = int_sorted[m-1:n-1]

    tmin = np.zeros((n,))
    tmin[:m] = int_sorted[0]
    tmin[m:n] = int_sorted[1:n-m+1]

    trez = np.array([v if v > np.power(10., -10) else np.power(10., -10) for v in (tadd - tmin)])
    vasr = 1/n*(np.sum(np.log(1.*n/(2*m)*trez)))+bi2_vas
    
    vasrnu = vasr-np.log(int_mean)

    return vasrnu


def calc_burst_index(intervals, bbh=DEFAULT_BBH, bbl=DEFAULT_BBL, brep=DEFAULT_REPEAT):
    bg1 = np.array(intervals <= bbh, dtype=int)
    bg2 = np.array(intervals <= bbl, dtype=int)

    bg123 = bg1 + bg2
    bgr = np.zeros((len(intervals),), dtype=int)

    for t in range(bg123.shape[0] - brep):
        if bg123[t] == 2:
            if np.sum(bg123[t:t+brep+1]) == 2*(brep+1):
                bgr[t:t+brep+1] = 2
            else:
                bgr[t] -= 1
        elif bg123[t] == 1:
            if np.sum(bg123[t:t+brep+1]) == 2*(brep+1) - 1:
                bgr[t:t+brep+1] = 2
            else:
                bgr[t] -= 1
        elif bg123[t] == 0:
            bgr[t] = 0

    bgrf = np.sum(bgr > 0)
    bgl2 = np.sum(intervals >= bbh)
    bil = bgrf/bgl2

    return bil


def calc_bi_two(spikes):
    int_1 = calc_intervals(spikes, step=1)
    int_2 = calc_intervals(spikes, step=2)
    e_1 = np.mean(int_1)

    return (2*np.var(int_1) - np.var(int_2))/(2*e_1*e_1)


def calc_freq_var(spikes):
    spikes = np.array(spikes)
    res = list()

    win_len = (spikes[-1] - spikes[0])/5
    step_size = win_len/5

    for s, e in [(s, s+win_len) for s in np.arange(np.floor(spikes[0]), np.ceil(spikes[-1] - win_len)+0.01, step_size)]:
        count = np.count_nonzero((spikes >= s)&(spikes <= e))
        res.append(1.*count/(e-s))

    return 100.*(max(res) - min(res))/max(res)


def calc_modalirity_burst(intervals, bound=DEFAULT_MODALIRITY_BOUND):
    below = np.sum(intervals < ms_to_sec(bound))
    above = np.sum(intervals >= ms_to_sec(bound))

    if above == 0:
        return 1.
    else:
        return below/above


def calc_pause_index(intervals, bound=DEFAULT_PAUSE_INDEX_BOUND):
    below = np.sum(intervals < ms_to_sec(bound))
    above = np.sum(intervals >= ms_to_sec(bound))

    if below == 0:
        return 1.
    else:
        return above/below


def calc_pause_ratio(intervals, bound=DEFAULT_PAUSE_RATIO_BOUND):
    above = np.sum(intervals >= ms_to_sec(bound)) + np.sum(intervals[intervals >= ms_to_sec(bound)])
    below = np.sum(intervals < ms_to_sec(bound)) + np.sum(intervals[intervals < ms_to_sec(bound)])

    if below == 0:
        return 1.
    else:
        return above/below


# what to do with low freq data? 
def calc_burst_behavior(intervals, spikes_count, bound=DEFAULT_BURST_BEHAVIOUR_BOUND):
    bound = ms_to_sec(bound)
    total_count = int(np.ceil(np.sum(intervals)/bound))
    above_count = 0

    partial_sum = np.cumsum(intervals)

    idx = 0
    n = len(partial_sum)
    for i in range(total_count):
        counter = 0
        upper_bound = (i+1)*bound
        while idx < n and partial_sum[idx] <= upper_bound:
            counter += 1
            idx += 1

        if counter >= spikes_count:
            above_count += 1


    return above_count/total_count


def calc_local_variance(isi):
    isi = np.asarray(isi)
    return 3.*np.sum(np.power(np.diff(isi)/(isi[:-1] + isi[1:]), 2))/(len(isi)-1)


def calc_burst_by_mean(intervals):
    return np.median(intervals)/np.mean(intervals)

    
def calc_skewness(intervals):
    return skew(intervals)


def calc_kurtosis(intervals):
    return kurtosis(intervals, fisher=False)

    
def get_type(med_mean, cv):
    if med_mean < 0.7:
        return 'burst'
    elif cv < 0.85:
        return 'tonic'
    else:
        return 'irregular'


def calc_isp(isi, hz_low, hz_high):
    f, Pxx_den = periodogram(isi, DEFAULT_FREQUENCY)
    idx_low = idx_of_nearest(f, hz_low)
    idx_high = idx_of_nearest(f, hz_high)
    
    return np.trapz(Pxx_den[idx_low:idx_high+1])
    

def calc_burst_rate(isi):
    isi = np.array(isi)
    return 1.*np.count_nonzero(isi > np.mean(isi))/len(isi)
