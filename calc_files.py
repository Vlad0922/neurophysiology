# -*- coding: utf-8 -*-

from __future__ import division

import sys
import os

import neo.io

from utility import *
from statistics import *
from oscore import *
 
ISP_RANGE = [(1,3), (3,8), (8,13), (13,30), (30,100)]
OSCORE_RANGE = [(3.,8.), (8.,12.), (12.,20.), (20.,30.), (30.,60.), (60., 90.)]

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
        Trial = sec_to_timestamps(data_filtered, DEFAULT_FREQUENCY).tolist()
        iTrialLength = Trial[-1]
        oscore = oscore_spikes(np.array([Trial]), iTrialLength, osc_l, osc_h, DEFAULT_FREQUENCY)
        df['oscore_{}_{}'.format(osc_l, osc_h)] = oscore
    
    return df

def main():
    dist_dir = sys.argv[1]
    dist_file = sys.argv[2] + '.xls'
      
    print 'main'
    for root, subdirs, files in os.walk(dist_dir):
        for full_name, f_name in [(root + '\\' + f_name, f_name) for f_name in files]:
            ext = full_name[-3:].lower()
            print full_name
            if ext == 'smr':
                r = neo.io.Spike2IO(filename=full_name)
                blks = r.read(cascade=True, lazy=False)
                for blk in blks:
                    for seg in blk.segments:
                        for st in seg.spiketrains:
                            spikes = np.array(st)
                            if len(spikes) > 50:
                                df = calc_stats(spikes, 'SMR neuron dummy', f_name)
                                write_to_excel(dist_file, 'all_results', df, ['doc_name', 'data_name'])   
            elif ext == 'nex':
                r = neo.io.NeuroExplorerIO(filename=full_name)
                blks = r.read(cascade=True, lazy=False)
                for blk in blks:
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
                                            spikes_filtered = spikes[np.where( (spikes >= s) & (spikes <= e))]
                                            df = calc_stats(spikes_filtered, f_name, st.name, int_name)
                                            write_to_excel(dist_file, 'all_results', df, ['doc_name', 'data_name', 'interval_name'])                                    
                            elif name_lower.startswith('allfile'):
                                spikes = np.array(st)
                                df = calc_stats(spikes, f_name, st.name,'allfile')
                                write_to_excel(dist_file, 'all_results', df, ['doc_name', 'data_name', 'interval_name']) 

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'Wrong usage!'
        print 'Usage: python calc_files.py <directory> <result_fname>'
        exit(1)

    main()
    print 'done'