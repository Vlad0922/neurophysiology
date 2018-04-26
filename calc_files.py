# -*- coding: utf-8 -*-
import sys
import os

import neo.io

# import nolds

import pandas as pd

from utility import *
from statistics import *
from oscore import *

import argparse

from detect_bursts import detect_with_logisi, find_burst_bunches, detect_with_cma, detect_with_ps


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


def calc_stats(spikes):
    data_filtered = spikes
    time_int = calc_intervals(spikes)

    if len(time_int) == 0:
        raise 'Empty filter result!'

    df = dict()
    df['burst_index'] = calc_burst_index(time_int, bbh=DEFAULT_BBH, bbl=DEFAULT_BBL, brep=DEFAULT_REPEAT)
    df['cv'] = calc_cv(time_int)
    df['nu'] = calc_nu(time_int)
    df['frequency_variance'] = calc_freq_var(data_filtered)
    df['modalirity_burst'] = calc_modalirity_burst(time_int)
    df['pause_index'] = calc_pause_index(time_int)
    df['pause_ratio'] = calc_pause_ratio(time_int)
    df['skewness'] = calc_skewness(time_int)
    df['kurtoisis'] = calc_kurtosis(time_int)
    df['AI'] = calc_burst_by_mean(time_int)
    df['type'] = get_type(df['AI'], df['cv'])
    df['isi_mean'] = np.mean(time_int)
    df['isi_median'] = np.median(time_int)
    df['isi_std'] = np.std(time_int)
    df['spike_count'] = len(data_filtered)
    df['filter_length'] = (data_filtered[-1] - data_filtered[0])
    df['bi_2'] = calc_bi_two(data_filtered)
    df['lv'] = calc_local_variance(time_int)
    df['firing_rate'] = 1.*len(data_filtered)/df['filter_length']
    df['burst_behaviour'] = calc_burst_behavior(time_int, int(np.ceil(df['firing_rate']/10)))
    df['burst_percent'] = calc_burst_percent(time_int)
    
    # df['dfa'] = nolds.dfa(time_int)
    # df['sampen'] = nolds.sampen(time_int)
    # df['corr_dim'] = nolds.corr_dim(time_int, 2)
    # df['hurst_rs'] = nolds.hurst_rs(time_int)

    for (osc_l, osc_h) in OSCORE_RANGE:
        trial = sec_to_timestamps(data_filtered, DEFAULT_FREQUENCY).tolist()
        trial_len = trial[~0]
        oscore = oscore_spikes(np.array([trial]), trial_len, osc_l, osc_h, DEFAULT_FREQUENCY)
        df['oscore_{}_{}'.format(osc_l, osc_h)] = oscore

    burst_args = dict_to_tuple({'skewness': 0.75, 'min_spike_count': 4, 'logisi_cutoff':0.1, 'si_threshold':3})
    burst_mask, burst_bunches, burst_lens = detect_with_vitek(spikes, burst_args)

    df['mean_burst_len'] = np.mean(burst_lens)
    df['ratio_burst_time'] = 1.*sum(burst_lens)/(spikes[~0] - spikes[0])
    df['burst_spike_percent'] = 1.*np.sum(burst_mask)/len(spikes)
    df['mean_spikes_in_burst'] = calc_mean_spikes_in_burst(burst_bunches)
    df['mean_isi_in_burst'] = calc_mean_isi_in_burst(burst_bunches)
    df['median_isi_in_burst'] = calc_median_isi_in_burst(burst_bunches)
    df['interburst_interval'] = calc_interburst_interval(burst_bunches)

    return df


def main(args):
    dist_dir = args.data_dir
    dist_file = '{}.xls'.format(args.dist_file)

    all_data = defaultdict(list)

    for root, subdirs, files in os.walk(dist_dir):

        for full_name, f_name in [(os.path.join(root, f_name), f_name) for f_name in files]:
            patient = full_name.split(os.sep)[3]
            ext = full_name[-3:].lower()

            if ext == 'smr' and False:
                print(full_name)
                for st in spiketrains_iterator(neo.io.Spike2IO(filename=full_name)):
                    spikes = np.array(st)
                    if len(spikes) > 40 and (spikes[~0] - spikes[0] > 5.):
                        df = calc_stats(spikes, f_name, 'SMR neuron dummy', 'interval dummy')
                        df = {key: str(val) for key, val in df.items()}
                        df['patient'] = patient
                        write_to_excel(dist_file, 'all_results', df, ['doc_name', 'data_name', 'interval_name'])
            elif ext == 'nex':
                print(full_name)
                r = neo.io.NeuroExplorerIO(filename=full_name)
                for blk in r.read(cascade=True, lazy=False):
                    for seg in blk.segments:
                        for st in seg.spiketrains:
                            name_lower = str(st.name.lower())
                            if name_lower.startswith('fon'):
                                spikes = np.array(st)
                                for interval in seg.epochs:
                                    int_name = interval.annotations['channel_name'].lower()
                                    if name_lower.startswith(int_name):
                                        for s, d in zip(interval.times, interval.durations):
                                            e = s + d
                                            spikes_filtered = spikes[np.where((spikes >= s) & (spikes <= e))]
                                            if len(spikes_filtered) > 40 and (spikes_filtered[~0] - spikes_filtered[0]) > 5.:
                                                df = calc_stats(spikes_filtered)
                                                df['patient'] = patient
                                                df['data_name'] = st.name
                                                df['doc_name'] = f_name
                                                df['interval_name'] = int_name

                                                for k in df:
                                                    all_data[k].append(df[k])
                            elif name_lower.startswith('allfile') or name_lower.startswith('nw'):
                                spikes = np.array(st)
                                if len(spikes) > 75 and (spikes[~0] - spikes[0] > 5.):                                    
                                    df = calc_stats(spikes)
                                    df['patient'] = patient
                                    df['patient'] = patient
                                    df['data_name'] = st.name
                                    df['doc_name'] = f_name
                                    df['interval_name'] = 'allfile'

                                    for k in df:
                                        all_data[k].append(df[k])

    all_data = pd.DataFrame(all_data)
    all_data.to_excel(dist_file, index=False)

    print('done')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Detect bursts in spike data')

    parser.add_argument('--data_dir', type=str, required=True,
                        help='Input directory')
    parser.add_argument('--dist_file', type=str, required=True,
                        help='File with results')
    parser.add_argument('--si_thresh', type=float, help='S parameter for PS method')

    args = parser.parse_args()

    main(args)