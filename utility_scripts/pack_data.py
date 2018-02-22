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

header_dtype = np.dtype([
    ('riff', (bytes, 5)),
    ('file_length', np.int32),
    ('WAVEfmt', (bytes, 5)),
    ('chunk_length', np.int32),
    ('format_tag', np.short),
    ('channels', np.short),
    ('SamplesPerSec', np.int32),
    ('AvgBytesPerSec', np.int32),
    ('blockAlign', np.short),
    ('bitsPerSec', np.short),
    ('data', (bytes, 4))
    ], align=0)


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


hz_map = {'ng': 20000,
          'pf': 1000,
          'snd': 10000,
          'a':1}


def find_hz(ext):
    for key, hz in hz_map.items():
        if key in ext.lower():
            return hz

def read_file(fname):
    with open(fname, 'rb') as f:
        filename, file_extension = os.path.splitext(fname)
        hz = find_hz(file_extension)

        if 'a' in file_extension.lower():
            f.seek(100)
            ftype = 'spikes'
            raw_data = np.fromfile(f, dtype=np.int32)
        else:
            ftype = 'signal'
            raw_data = np.fromfile(f, dtype=np.int16)

    return raw_data, ftype, hz


def write_to_smr(data_dict, fname):
    max_hz = max([hz for (_, _, hz) in data_dict.values()])

    fullpath_smr = fname
    print(fname)

    fhand = ceds_handler.S64Create(fullpath_smr, 32, 0)
    ceds_handler.S64SetTimeBase(fhand, ctypes.c_double(1./max_hz)) 

    for idx, key in enumerate(data_dict.keys(), start=1):
        arr, ftype, hz = data_dict[key]
        
        div = max_hz//hz

        ceds_handler.S64SetIdealRate(fhand, ctypes.c_int(idx), ctypes.c_double(max_hz))

        if ftype == 'signal':
            arr = [ctypes.c_short(v) for v in arr]
            ceds_handler.S64SetWaveChan(ctypes.c_int(fhand), ctypes.c_int(idx), ctypes.c_longlong(div), ctypes.c_int(1), ctypes.c_double(hz))
            ceds_handler.S64SetChanScale(fhand, idx, ctypes.c_double(6553.6))

            res_write = ceds_handler.S64WriteWaveS(ctypes.c_int(fhand), ctypes.c_int(idx), (ctypes.c_short * len(arr))(*arr), 
                                 ctypes.c_int(len(arr)), ctypes.c_longlong(0))

            ceds_handler.S64SetChanUnits(fhand, idx, 'mV')
        else:
            arr = [ctypes.c_longlong(v) for v in arr]
            set_ev_res = ceds_handler.S64SetEventChan(ctypes.c_int(fhand), ctypes.c_int(idx), ctypes.c_double(max_hz), ctypes.c_int(2))
            ev_res = ceds_handler.S64WriteEvents(ctypes.c_int(fhand), ctypes.c_int(idx), (ctypes.c_longlong * len(arr))(*arr), ctypes.c_int(len(arr)))

        ceds_handler.S64SetChanTitle(fhand, idx, key)
        
        
        

    ceds_handler.S64Close(fhand)


def main(args):
    dirname = args.data_dir
    data_dct = dict()
    for fname in [f for f in os.listdir(dirname) if os.path.isfile(os.path.join(dirname, f))]:
        full_path = os.path.join(dirname, fname)
        data_dct[fname.replace('.', '_')] = read_file(full_path)


    print(data_dct.keys())

    write_to_smr(data_dct, dirname + '.smr')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert AlphaOmega .map files to .txt/.smr for spike2')

    parser.add_argument('--data_dir', type=str, required=True,
                        help='Directory with data files')

    args = parser.parse_args()
    main(args)