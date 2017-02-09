from __future__ import division


import numpy as np
from scipy.special import psi


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

def sec_to_ms(sec):
    return sec*1000


def ms_to_sec(ms):
    return ms/1000


def calc_intervals(spike_data, filter_by=DEFAULT_FILT_TYPE, low=DEFAULT_FILT_LOWER, high=DEFAULT_FILT_UPPER):
    if low == -1 or high == -1:
        res = np.ediff1d(spike_data)
    elif filter_by == FILTER_TIME:
        spike_data = spike_data[(spike_data >= low) & (spike_data <= high)]
        res = np.ediff1d(spike_data)
    elif filter_by == FILTER_NUMBER:
        res = np.ediff1d(spike_data)[low:high+1]
    else:
        raise Exception('Wrong filter type!')

    return res


def calc_cv(intervals):
    n = len(intervals)
    int_mean = np.mean(intervals)
    int_var = np.var(intervals)
    int_cv = np.sqrt(int_var)/int_mean

    return int_cv


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

    trez = tadd - tmin
    vasr = 1/n*(np.sum(np.log(n/(2*m)*trez)))+bi2_vas
    vasrnu = vasr-np.log(int_mean)

    return vasrnu


def calc_bi(intervals, bbh=DEFAULT_BBH, bbl=DEFAULT_BBL, brep=DEFAULT_REPEAT):
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


def calc_moment(data, moment):
    m = np.mean(data)
    def sum_moment(v):
        return np.power(v - m, moment)

    t = np.sum(np.apply_along_axis(sum_moment, 0, data))
    return t/((len(data) - 1)*np.power(np.std(data), moment))  


def calc_freq_var(intervals):
    return 100*(max(intervals) - min(intervals))/max(intervals)


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
def calc_burst_behavior(intervals, bound=DEFAULT_BURST_BEHAVIOUR_BOUND):
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

        if counter >= 5:
            above_count += 1


    return above_count/total_count


def calc_skewness(intervals):
    return calc_moment(intervals, 3)


def calc_kurtosis(intervals):
    return calc_moment(intervals, 4)
