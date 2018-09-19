# -*- coding: utf-8 -*-
import sys
import os
import glob

import neo.io

import pandas as pd

from statistics import *
from oscore import *

import argparse

from detect_bursts import detect_with_ps, detect_with_vitek

from collections import namedtuple, defaultdict


ISP_RANGE = [(1, 3), (3, 8), (8, 13), (13, 30), (30, 100)]


def dict_to_tuple(dictionary):
    return namedtuple('GenericDict', dictionary.keys())(**dictionary)


def spiketrains_iterator(handler):
    for blk in handler.read(cascade=True, lazy=False):
        for seg in blk.segments:
            for st in seg.spiketrains:
                yield st

    raise StopIteration


def safe_mean(arr):
    if len(arr) == 0:
        return 0.
    else:
        return np.mean(arr)


def calc_stats(spikes, args):
    if args.burst_algo == 'PS':
        burst_method = detect_with_ps
    elif args.burst_algo == 'vitek':
        burst_method = detect_with_vitek

    data_filtered = spikes
    isi = calc_intervals(spikes)

    if len(isi) == 0:
        raise RuntimeError('Empty filter result!')

    df = dict()
    df['burst_index'] = calc_burst_index(isi, bbh=DEFAULT_BBH, bbl=DEFAULT_BBL, brep=DEFAULT_REPEAT)
    df['cv'] = calc_cv(isi)
    df['diff_entropy (Nu)'] = calc_nu(isi)
    df['frequency_variance'] = calc_freq_var(data_filtered)
    df['modalirity_burst'] = calc_modalirity_burst(isi)
    df['pause_index'] = calc_pause_index(isi)
    df['pause_ratio'] = calc_pause_ratio(isi)
    df['skewness'] = calc_skewness(isi)
    df['kurtoisis'] = calc_kurtosis(isi)
    df['AI'] = calc_burst_by_mean(isi)
    df['type'] = get_type(df['AI'], df['cv'])
    df['isi_mean'] = np.mean(isi)
    df['isi_median'] = np.median(isi)
    df['isi_std'] = np.std(isi)
    df['spike_count'] = len(data_filtered)
    df['filter_length'] = (data_filtered[~0] - data_filtered[0])
    df['bi_2'] = calc_bi_two(data_filtered)
    df['local_variance'] = calc_local_variance(isi)
    df['firing_rate'] = 1.*len(data_filtered)/df['filter_length']
    df['burst_behaviour'] = calc_burst_behavior(isi, int(np.ceil(df['firing_rate']/10)))
    df['burst_percent'] = calc_burst_percent(isi)
    df['ISI_larger_mean'] = 1.*sum(isi >= df['isi_mean'])/len(isi)

    for (osc_l, osc_h) in OSCORE_RANGE:
        trial = sec_to_timestamps(data_filtered, DEFAULT_FREQUENCY).tolist()
        trial_len = trial[~0]
        oscore = oscore_spikes(np.array([trial]), trial_len, osc_l, osc_h, DEFAULT_FREQUENCY)
        df['oscore_{}_{}'.format(osc_l, osc_h)] = oscore

    burst_args = dict_to_tuple({'skewness': 0.75, 'min_spike_count': 4, 'logisi_cutoff':0.1, 'si_threshold':3})
    burst_mask, burst_bunches, burst_lens = burst_method(spikes, burst_args)

    df['burst_rate'] = 1.*len(burst_lens)/(spikes[~0] - spikes[0])
    df['mean_burst_len'] = np.mean(burst_lens)
    df['ratio_burst_time'] = 1.*sum(burst_lens)/(spikes[~0] - spikes[0])
    df['burst_spike_percent'] = 1.*np.sum(burst_mask)/len(spikes)
    df['mean_spikes_in_burst'] = calc_mean_spikes_in_burst(burst_bunches)
    df['mean_isi_in_burst'] = calc_mean_isi_in_burst(burst_bunches)
    df['median_isi_in_burst'] = calc_median_isi_in_burst(burst_bunches)
    df['interburst_interval'] = calc_interburst_interval(burst_bunches)

    return df


def merge_st(st_list):
    print(st_list)
    offset = st_list[0][~0];

    for i in range(1, len(st_list)):
        st_list[i] = st_list[i][1:] - st_list[i][0] + offset
        offset = st_list[i][~0]

    return np.concatenate(st_list)


def main(args):
    print(args.all)
    dist_dir = args.data_dir
    dist_file = '{}.xls'.format(args.dist_file)

    all_data = defaultdict(list)

    for root, subdirs, files in os.walk(dist_dir):
        for full_name, f_name in [(os.path.join(root, f_name), f_name) for f_name in files]:
            patient = full_name.split(os.sep)[~1]

            ext = full_name.split('.')[~0].lower()

            if not('nex' in ext):
                continue

            print(full_name)
            r = neo.io.NeuroExplorerIO(filename=full_name)
            for blk in r.read(cascade=True, lazy=False):
                for seg in blk.segments:
                    for st in seg.spiketrains:
                        name_lower = str(st.name.lower())
                        if name_lower.startswith('fon') or args.all:
                            spikes = np.array(st)
                            for interval in seg.epochs:
                                int_name = interval.annotations['channel_name'].lower()
                                if name_lower.startswith(int_name) or args.all:
                                    print('\t {:15} \t {:15}'.format(name_lower, int_name))
                                    int_spikes = list()
                                    interval_names = list()
                                    for s, d in zip(interval.times, interval.durations):
                                        e = s + d
                                        spikes_filtered = spikes[np.where((spikes >= s) & (spikes <= e))]
                                        print(s, e, spikes[0], spikes[~0])

                                        print(len(spikes_filtered), spikes_filtered[~0] - spikes_filtered[0])
                                    
                                        if len(spikes_filtered) > 50:
                                            int_spikes.append(spikes_filtered)
                                            interval_names.append(int_name)

                                    spikes_merged = merge_st(int_spikes) # часть алгоритмов (oscore) требует спайков, а не иси, поэтому смержим ST
                                    df = calc_stats(spikes_merged, args)
                                    df['patient'] = patient
                                    df['data_name'] = st.name
                                    df['doc_name'] = f_name
                                    df['interval_name'] = ','.join(interval_names)

                                    for k in df:
                                        all_data[k].append(df[k])
                        elif name_lower.startswith('allfile'):
                            spikes = np.array(st)
                            if len(spikes) > 50 and (spikes[~0] - spikes[0] > 5.):
                                df = calc_stats(spikes, args)
                                df['patient'] = patient
                                df['data_name'] = st.name
                                df['doc_name'] = f_name
                                df['interval_name'] = 'allfile'

                                for k in df:
                                    all_data[k].append(df[k])

    all_data = pd.DataFrame(all_data)
    all_data = all_data[['patient', 'doc_name', 'data_name', 'interval_name', 'filter_length', 'spike_count', 'type', 'firing_rate', 'cv', 'AI', 'frequency_variance',  'isi_mean', 'isi_median',
                          'isi_std', 'skewness', 'kurtoisis', 'local_variance',  'diff_entropy (Nu)', 'ISI_larger_mean', 'burst_index', 'burst_spike_percent', 'ratio_burst_time', 'burst_rate',  
                          'interburst_interval',  'mean_burst_len',  'mean_isi_in_burst', 'mean_spikes_in_burst', 'median_isi_in_burst', 'pause_index',
                          'oscore_3.0_8.0',  'oscore_8.0_12.0', 'oscore_12.0_20.0', 'oscore_20.0_30.0', 'oscore_30.0_60.0', 'oscore_60.0_90.0']];

    all_data.to_excel(dist_file, index=False)

    print('done')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Detect bursts in spike data')

    parser.add_argument('--data_dir',   type=str, required=True, help='Input directory')
    parser.add_argument('--dist_file',  type=str, required=True, help='File with results')
    parser.add_argument('--si_thresh',  type=float, help='S parameter for PS method')
    parser.add_argument('--bin_func',   type=str, default='discharge')
    parser.add_argument('--burst_algo', type=str, default='PS')
    parser.add_argument('--all', type=bool, default=False)

    args = parser.parse_args()

    main(args)
