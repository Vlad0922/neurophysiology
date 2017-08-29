# -*- coding: utf-8 -*-

from __future__ import division

import sys
import os

import neo.io

from utility import *
from statistics import *
from oscore import *

import argparse

from detect_bursts import detect_with_logisi, find_burst_bunches



ISP_RANGE = [(1, 3), (3, 8), (8, 13), (13, 30), (30, 100)]


def dict_to_tuple(dictionary):
    return namedtuple('GenericDict', dictionary.keys())(**dictionary)


def spiketrains_iterator(handler):
    # try:
    for blk in handler.read(cascade=True, lazy=False):
        for seg in blk.segments:
            for st in seg.spiketrains:
                yield st
    # except:
    #     print 'Something wrong with file:'.format(handler.filename)
    #     raise StopIteration

    raise StopIteration


def calc_stats(spikes, fname, neuron_name, interval_name):
    data_filtered = spikes
    time_int = calc_intervals(spikes)

    if len(time_int) == 0:
        raise 'Empty filter result!'

    df = dict()
    df['data_name'] = neuron_name
    df['burst_index'] = calc_burst_index(time_int, bbh=DEFAULT_BBH, bbl=DEFAULT_BBL, brep=DEFAULT_REPEAT)
    df['cv'] = calc_cv(time_int)
    df['nu'] = calc_nu(time_int)
    df['frequency_variance'] = calc_freq_var(data_filtered)
    df['modalirity_burst'] = calc_modalirity_burst(time_int)
    df['pause_index'] = calc_pause_index(time_int)
    df['pause_ratio'] = calc_pause_ratio(time_int)
    df['skewness'] = calc_skewness(time_int)
    df['kurtoisis'] = calc_kurtosis(time_int)
    df['burst_mean'] = calc_burst_by_mean(time_int)
    df['type'] = get_type(df['burst_mean'], df['cv'])
    df['doc_name'] = fname
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
    df['interval_name'] = interval_name

    for (osc_l, osc_h) in OSCORE_RANGE:
        trial = sec_to_timestamps(data_filtered, DEFAULT_FREQUENCY).tolist()
        trial_len = trial[~0]
        oscore = oscore_spikes(np.array([trial]), trial_len, osc_l, osc_h, DEFAULT_FREQUENCY)
        df['oscore_{}_{}'.format(osc_l, osc_h)] = oscore

    burst_args = dict_to_tuple({'skewness': 0.75, 'min_spike_count': 4, 'logisi_cutoff':0.1})
    burst_mask, burst_bunches = detect_with_logisi(spikes, burst_args)

    df['burst_spike_percent'] = 1.*np.sum(burst_mask)/len(spikes)
    df['mean_spikes_in_burst'] = calc_mean_spikes_in_burst(burst_bunches)
    df['mean_isi_in_burst'] = calc_mean_isi_in_burst(burst_bunches)
    df['median_isi_in_burst'] = calc_median_isi_in_burst(burst_bunches)
    df['interburst_interval'] = calc_interburst_interval(burst_bunches)

    bursts_oscore = calc_oscore_for_bursts(burst_bunches)
    df.update(bursts_oscore)

    return df


def main(args):
    dist_dir = args.data_dir
    dist_file = '{}.xls'.format(args.dist_file)

    for root, subdirs, files in os.walk(dist_dir):
        for full_name, f_name in [(os.path.join(root, f_name), f_name) for f_name in files]:
            patient = full_name.split(os.sep)[3]
            ext = full_name[-3:].lower()

            if ext == 'smr':
                print full_name
                for st in spiketrains_iterator(neo.io.Spike2IO(filename=full_name)):
                    spikes = np.array(st)
                    if len(spikes) > 40 and (spikes[~0] - spikes[0] > 5.):
                        df = calc_stats(spikes, f_name, 'SMR neuron dummy', 'interval dummy')
                        df = {key: str(val) for key, val in df.items()}
                        df['patient'] = patient
                        write_to_excel(dist_file, 'all_results', df, ['doc_name', 'data_name', 'interval_name'])
            elif ext == 'nex':
                print full_name
                r = neo.io.NeuroExplorerIO(filename=full_name)
                for blk in r.read(cascade=True, lazy=False):
                    for seg in blk.segments:
                        for st in seg.spiketrains:
                            name_lower = st.name.lower()
                            if name_lower.startswith('fon'):
                                spikes = np.array(st)
                                for interval in seg.epochs:
                                    int_name = interval.annotations['channel_name'].lower()
                                    if name_lower.startswith(int_name):
                                        for s, d in zip(interval.times, interval.durations):
                                            e = s + d
                                            spikes_filtered = spikes[np.where((spikes >= s) & (spikes <= e))]
                                            if len(spikes_filtered) > 40 and (spikes_filtered[~0] - spikes_filtered[0]) > 5.:
                                                df = calc_stats(spikes_filtered, f_name, st.name, int_name)
                                                df['patient'] = patient
                                                write_to_excel(dist_file, 'all_results', df, ['doc_name', 'data_name', 'interval_name'])
                            elif name_lower.startswith('allfile'):
                                spikes = np.array(st)
                                if len(spikes) > 40 (spikes[~0] - spikes[0] > 5.):                                    
                                    df = calc_stats(spikes, f_name, st.name, 'allfile')
                                    df['patient'] = patient
                                    write_to_excel(dist_file, 'all_results', df, ['doc_name', 'data_name', 'interval_name'])

    print 'done'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Detect bursts in spike data')

    parser.add_argument('--data_dir', type=str, required=True,
                        help='Input directory')
    parser.add_argument('--dist_file', type=str, required=True,
                        help='File with results')

    args = parser.parse_args()

    main(args)