# -*- coding: utf-8 -*-

import sys
import os
import glob
import ctypes
import argparse
import platform

import math

from collections import defaultdict

import neo.io

import numpy as np


BYTES_IN_KB = 1024
OK_SIGNALS_COUNT = 3


def get_windows_version():
    if platform.architecture()[0] == '64bit':
        return 'x64'
    else:
        return 'x86'


def get_bin_path():
    curr_dir = os.path.split(os.getcwd())[~0]
    if curr_dir == 'neurophysiology':
        return os.path.join('utility_scripts', 'ceds_bin', WIN_VER)
    else:
        return os.path.join('ceds_bin', WIN_VER)


WIN_VER = get_windows_version()
CEDS_BIN_PATH = get_bin_path()


ceds_path = os.path.join(CEDS_BIN_PATH, 'ceds64int.dll')
son_path =  os.path.join(CEDS_BIN_PATH, 'son64.dll')
msvcp_path =  os.path.join(CEDS_BIN_PATH, 'msvcp110.dll')
msvcr_path =  os.path.join(CEDS_BIN_PATH, 'msvcr110.dll')
load_order = [msvcr_path, msvcp_path, son_path]


def load_ceds_lib():
    for l in load_order:
        ctypes.WinDLL(l)

    return ctypes.WinDLL(ceds_path)
ceds_handler = load_ceds_lib()


def read_alphaomega(fname):
    data_dict = dict()

    r = neo.io.AlphaOmegaIO(filename=fname)
    for blk in r.read(cascade=True, lazy=False):
        for seg in blk.segments:
            lfp_signals = sorted(sg.name for sg in seg.analogsignals if sg.name.lower().startswith('lfp'))[:OK_SIGNALS_COUNT]
            spk_signals = sorted(sg.name for sg in seg.analogsignals if sg.name.lower().startswith('spk'))[:OK_SIGNALS_COUNT]
            for sg in seg.analogsignals:
                if sg.name in lfp_signals or sg.name in spk_signals:
                    data_dict[sg.name] = [v[0] for v in np.array(sg)]

    return data_dict


def rolling_window(a, window):
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def _interpolate_window(w, ratio):
    res = np.zeros(ratio)
    res[0], res[~0] = w

    for i in range(1, ratio):
        res[i] = w[0] + 1.*(w[1] - w[0])*i/(ratio+1)

    return res 


def inteprolate_signal(signal, ratio):
    windows = rolling_window(signal, 2)
    return np.hstack([_interpolate_window(w, ratio) for w in windows])


def prepare_keys(data_dict):
    orig_keys = data_dict.keys()
    for k in orig_keys:
        data_dict[k.replace(' ', '_')] = data_dict.pop(k)

    return data_dict


def write_to_smr(data_dict, root, fname, spk_freq, lfp_freq):
    div = int(spk_freq/lfp_freq)

    preproc_path = os.path.join(root, 'converted_to_smr')
    if not os.path.exists(preproc_path):
        os.makedirs(preproc_path)

    fname_no_ext = os.path.splitext(os.path.basename(fname))[0]
    fname_smr = '{}.smr'.format(fname_no_ext)
    fullpath_smr = os.path.join(preproc_path, fname_smr)   

    fhand = ceds_handler.S64Create(fullpath_smr, 32, 0)
    ceds_handler.S64SetTimeBase(fhand, ctypes.c_double(1./spk_freq)) 

    for idx, key in enumerate(data_dict.keys()):
        idx += 1

        arr = [ctypes.c_float(v) for v in data_dict[key]]

        ceds_handler.S64SetIdealRate(fhand, ctypes.c_int(idx), ctypes.c_double(spk_freq))

        if key.lower().startswith('lfp'):
            ceds_handler.S64SetWaveChan(ctypes.c_int(fhand), ctypes.c_int(idx), ctypes.c_longlong(div), ctypes.c_int(9), ctypes.c_double(lfp_freq))
            ceds_handler.S64SetChanTitle(fhand, idx, key)
        else:
            ceds_handler.S64SetWaveChan(ctypes.c_int(fhand), ctypes.c_int(idx), ctypes.c_longlong(1), ctypes.c_int(9), ctypes.c_double(spk_freq))
            ceds_handler.S64SetChanTitle(fhand, idx, key)

        res_write = ceds_handler.S64WriteWaveF(ctypes.c_int(fhand), ctypes.c_int(idx), (ctypes.c_float * len(arr))(*arr), 
                                         ctypes.c_int(len(arr)), ctypes.c_longlong(0))

        
        ceds_handler.S64SetChanUnits(fhand, idx, 'mV')

    ceds_handler.S64Close(fhand)


def get_files_by_trials(dat_files):
    files_by_trial = defaultdict(lambda: list())

    for full_path in dat_files:
        fname = os.path.basename(full_path)
        num = fname[:2]
        files_by_trial[num].append(full_path)

    return files_by_trial


def read_binary(files, cut, cut_size, freq):
    data_dict = dict()
    for full_path in files:
        fname = os.path.splitext(os.path.basename(full_path))[0]
        fname = fname[fname.find('_')+1:]

        data = np.fromfile(full_path, dtype="<u2")

        if cut:
            data = data[int(cut_size*freq):]

        data_dict[fname] = data

    return data_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert AlphaOmega .map files to .txt/.smr for spike2')

    parser.add_argument('--data_dir', type=str, required=True,
                        help='Directory with data files')
    parser.add_argument('--spk_freq', type=float, default=24340.7688141,
                        help='Spike channel frequency')
    parser.add_argument('--lfp_freq', type=float, default=760.648965836,
                        help='LFP channel frequency')
    parser.add_argument('--cut', action='store_true', default=False,
                        help='Cut strange artifact at the beginning?')
    parser.add_argument('--cut_size', type=float, help='Cut size')

    args = parser.parse_args()

    if args.cut and args.cut_size is None:
        parser.error("--cut requires --cut_size")

    input_dir = args.data_dir
    lfp_freq = args.lfp_freq
    spk_freq = args.spk_freq
    cut = args.cut
    cut_size = args.cut_size

    walk_res = os.walk(input_dir)

    for root, _, _ in walk_res:
        map_files = glob.glob(os.path.join(root, '*.map'))
        dat_files = [fname for fname in glob.glob(os.path.join(root, '*.dat')) if os.path.getsize(fname) > 2*BYTES_IN_KB]

        for fname in map_files:
            print 'Working with: {}'.format(fname)

            data_dict = read_alphaomega(fname)
            data_dict = prepare_keys(data_dict)

            write_to_smr(data_dict, root, fname, spk_freq, lfp_freq)

        files_by_trial = get_files_by_trials(dat_files)
        for num, files in files_by_trial.items():
            print 'Working with: {}'.format(root)

            data_dict = read_binary(files, cut, cut_size, spk_freq)
            data_dict = prepare_keys(data_dict)

            write_to_smr(data_dict, root, str(num), spk_freq, lfp_freq)

    print 'done'