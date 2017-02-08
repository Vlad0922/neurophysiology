from __future__ import division


import numpy as np
from scipy.special import psi


FILTER_TIME = 1
FILTER_NUMBER = 0

def sec_to_ms(sec):
    return sec*1000


def ms_to_sec(ms):
    return ms/1000


def calc_intervals(spike_data, filter_by, low, high):
    if filter_by == FILTER_TIME:
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


def calc_nu(intervals, wsize):
    n = len(intervals)
    m = wsize
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


def calc_bi(intervals, bbh, bbl, brep):
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
    m = np.mean(intervals)
    def sum_moment(v, m):
        return np.power(v - m, moment)

    t = np.sum(np.apply_along_axis(sum_moment, 0, intervals))
    return t/((len(intervals) - 1)*np.power(np.std(intervals), moment))  


def calc_freq_var(intervals):
    return 100*(max(intervals) - min(intervals))/max(intervals)


def calc_modalirity_burst(intervals):
    return len(intervals < ms_to_sec(10))/len(intervals >= ms_to_sec(10))

def calc_pause_index(intervals):
    return len(intervals > ms_to_sec(50))/len(intervals <= ms_to_sec(50))


def calc_pause_ratio(intervals):
    return (len(intervals > ms_to_sec(50)) + np.sum(intervals[intervals > ms_to_sec(50)]))/(len(intervals <= ms_to_sec(50)) + np.sum(intervals[intervals <= ms_to_sec(50)]))


def calc_burst_behavior(intervals):
    bound = ms_to_sec(100)
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


def calc_kurtosis(intevals):
    return calc_moment(intervals, 4)
