from __future__ import division

import imp
import sys
import re

statistics = imp.load_source('statistics', '../../../common/statistics.py')
from statistics import *

import numpy as np
import pandas as pd

PARAMS = dict()

def frame_to_time(spikes, freq):
    return np.array([s/freq for s in spikes])


def read_freq(info_line):
    frame_dat = info_line[info_line.find('(')+1:info_line.find(')')] 

    frames, _, _, time, measure = frame_dat.split()

    frames = int(frames)
    time = float(time)

    if measure == 'sec':
        time *= 1000

    return int(1000/time/frames)

def read_gaya_file(filename):
    mode = PARAMS['mode']
    with open(filename) as f:
        print 'processing ', filename, '...'
        if mode == 'temp':
            temp = int(re.findall(r"[\w']+", filename.replace('_', ','))[4])

        file_info = f.readline()
        freq = read_freq(file_info)
        param_coeff = 1000/freq

        f.readline() # pass empty line
        headers = f.readline()

        df = pd.DataFrame(columns=['name', 'x', 'y', 'burst_index', 'cv', 'nu', 'freq_variance', 'modalirity_burst',
                                    'pause_index', 'pause_ratio', 'burst_behavior', 'skewness', 'kurtosis'])

        for line in f.readlines():
            metadata, spikes = line.split(':')
            cell, x, y, n_spikes = metadata.split()
            n_spikes = int(n_spikes)
            x = float(x)
            y = float(y)

            if n_spikes >= DEFAULT_WIN_SIZE:
                spikes = frame_to_time(map(int, spikes.split()), freq)
                time_int = calc_intervals(spikes)

                bi = calc_burst_index(time_int, bbh=param_coeff*DEFAULT_BBH, bbl=param_coeff*DEFAULT_BBL)
                cv = calc_cv(time_int)
                nu = calc_nu(time_int)
                freq_v = calc_freq_var(time_int)
                mod_burst = calc_modalirity_burst(time_int, param_coeff*DEFAULT_MODALIRITY_BOUND)
                pause_ind = calc_pause_index(time_int, param_coeff*DEFAULT_PAUSE_INDEX_BOUND)
                pause_rat = calc_pause_ratio(time_int, param_coeff*DEFAULT_PAUSE_RATIO_BOUND)
                burst_beh = calc_burst_behavior(time_int, param_coeff*DEFAULT_BURST_BEHAVIOUR_BOUND)
                skew = calc_skewness(time_int)
                kurt = calc_kurtosis(time_int)

                df = df.append({'name':cell, 'x':x, 'y':y, 'burst_index':bi, 'cv':cv, 'nu':nu, 'freq_variance':freq_v, 'modalirity_burst':mod_burst,
                                'pause_index':pause_ind, 'pause_ratio':pause_rat, 'burst_behavior':burst_beh, 'skewness':skew, 'kurtosis':kurt}, ignore_index=True)

        if mode == 'temp':
            df['temp'] = temp

    return df


def main():
    files_list = sys.argv[2:-1]

    df = pd.DataFrame(columns=['name', 'x', 'y', 'burst_index', 'cv', 'nu', 'freq_variance', 'modalirity_burst',
                            'pause_index', 'pause_ratio', 'burst_behavior', 'skewness', 'kurtosis'])
    for f_name in files_list:
        file_data = read_gaya_file(f_name)
        df = pd.concat([df, file_data])

    df.to_csv(sys.argv[-1], index=False)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'wrong usage! Input file required\n'
        exit(1)

    PARAMS['mode'] = sys.argv[1]

    main()
