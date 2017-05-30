# -*- coding: utf-8 -*-

from __future__ import division

import itertools

from collections import defaultdict

import sys
import os

import neo.io

from statistics import *

BIN_SIZE = 0.2
WINDOW_SIZE = 1.

def parse_interval_file(fname):
    res = defaultdict(lambda: list())
    with open(fname) as f:
        for line in f.readlines():
            key, token, i_start, i_end = line.rstrip().split()
            res[key].append((int(i_start), int(i_end)))

    return res

def calc_stats(spikes):
    return (np.mean(spikes), np.median(spikes), len(spikes)/(spikes[-1] - spikes[0]))


def merge_intervals(intervals):
    res = list()
    curr_start, curr_end = intervals[0]
    for i_start, i_end in intervals[1:]:
        if np.abs(i_start - curr_end) > 0.005:
            res.append((curr_start, curr_end))
            curr_start, curr_end = i_start, i_end
        else:
            curr_end = i_end

    res.append((curr_start, curr_end))

    return res



def filter_len(intervals):
    return [inter for inter in intervals if (inter[1] - inter[0]) > 2.]

def main():
    data_file = sys.argv[1]
    interval_file = sys.argv[2]

    intervals = parse_interval_file(interval_file)

    r = neo.io.NeuroExplorerIO(filename=data_file)
    blks = r.read(cascade=True, lazy=False)
    for blk in blks:
        for seg in blk.segments:
            for st in seg.spiketrains:
                spikes = np.array(st)

                st_intervals = [0] + list(itertools.chain.from_iterable(intervals[st.name])) + [spikes[-1]]

                it = iter(st_intervals)
                st_intervals = zip(it, it)

                st_intervals = filter_len(st_intervals)

                res = list()
                for i_start, i_end, in st_intervals:
                    int_spikes = spikes[np.where((spikes >= i_start) & (spikes < i_end))]
                    for w_start, w_end in [(s, s+WINDOW_SIZE) for s in np.arange(i_start, i_end - WINDOW_SIZE + 0.01)]:
                        w_spikes = int_spikes[np.where((int_spikes >= w_start)&(int_spikes < w_end))]
                        if len(w_spikes) == 0:
                            break
                        w_params = list()
                        for b_start, b_end in [(b, b+BIN_SIZE) for b in np.arange(w_start, w_end, BIN_SIZE)]:
                            bin_spikes = w_spikes[np.where((w_spikes >= b_start)&(w_spikes < b_end))]
                            if len(bin_spikes) == 0:
                                break
                            w_params.append(calc_stats(bin_spikes))

                        if len(w_params) < 4:
                            break

                        params_mean = np.mean(w_params)
                        std = np.std(np.var(w_params, axis=0))

                        avg = params_mean/std
                        if np.mean(avg) < 0.33:
                            res.append((w_start, w_end))

                if len(res) == 0:
                    continue

                res = merge_intervals(res)
                

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print 'wrong usage!'
        print 'usage find_background.py <data_file> <interval_file>'
        exit(1)

    main()