# -*- coding: utf-8 -*-

import sys
import os
import glob

import neo.io

from collections import defaultdict

import numpy as np

import argparse

try:
    import matlab.engine

    eng = matlab.engine.start_matlab()
    eng.addpath('utility_scripts')

    USE_MATLAB = True
except:
    import pandas as pd

    print 'Cannot import or start MATLAB engine, will import to .txt files...'

    USE_MATLAB = False


BYTES_IN_KB = 1024
UV_IN_VOLTS = 0.001
OK_SIGNALS_COUNT = 3
KNOWN_TYPES = ('*.map', '*.dat')


def load_ceds_lib(lib_path):
    eng.addpath(lib_path, nargout=0)
    eng.CEDS64LoadLib(lib_path, nargout=0)


def read_alphaomega(fname):
    data_dict = dict()

    r = neo.io.AlphaOmegaIO(filename=fname)
    for blk in r.read(cascade=True, lazy=False):
        for seg in blk.segments:
            lfp_signals = sorted(sg.name for sg in seg.analogsignals if sg.name.lower().startswith('lfp'))[:OK_SIGNALS_COUNT]
            spk_signals = sorted(sg.name for sg in seg.analogsignals if sg.name.lower().startswith('spk'))[:OK_SIGNALS_COUNT]
            for sg in seg.analogsignals:
                if sg.name in lfp_signals or sg.name in spk_signals:
                    data_dict[sg.name] = [float(v[0]) for v in np.array(sg)]

    return data_dict


def prepare_keys(data_dict):
    orig_keys = data_dict.keys()
    for k in orig_keys:
        data_dict[k.replace(' ', '_')] = data_dict.pop(k)

    return data_dict


def write_with_matlab(data_dict, root, fname, spk_freq, lfp_freq):
    preproc_path = os.path.join(root, 'converted_to_smr')
    if not os.path.exists(preproc_path):
        os.makedirs(preproc_path)

    for k, v in data_dict.items():
        data_dict[k] = matlab.double(list(v))

    fname_no_ext = os.path.splitext(os.path.basename(fname))[0]
    fname_smr = '{}.smr'.format(fname_no_ext)
    fullpath_smr = os.path.join(preproc_path, fname_smr)

    eng.dict_to_smr(data_dict, fullpath_smr, spk_freq, lfp_freq, nargout=0)


def write_to_txt(data_dict, root, fname):
    preproc_path = os.path.join(root, 'converted_to_txt')
    if not os.path.exists(preproc_path):
        os.makedirs(preproc_path)

    df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in data_dict.iteritems()]))
    df = df[sorted(df.columns.values, reverse=True)]

    fname_no_ext = os.path.splitext(os.path.basename(fname))[0]
    fname_txt = '{}.txt'.format(fname_no_ext)
    fullpath_txt = os.path.join(preproc_path, fname_txt)

    df.to_csv(fullpath_txt, index=False, sep='\t')


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

        data = np.fromfile(full_path, dtype="<u2")*UV_IN_VOLTS

        if cut:
            data = data[int(cut_size*freq):]

        data_dict[fname] = data

    return data_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert AlphaOmega .map files to .txt/.smr for spike2')

    parser.add_argument('--input_dir', type=str, required=True,
                        help='Directory with data files')
    parser.add_argument('--spk_freq', type=float, required=True,
                        help='Spike channel frequency')
    parser.add_argument('--lfp_freq', type=float, required=True,
                        help='LFP channel frequency')
    parser.add_argument('--ceds_lib', type=str, default=r'C:\CEDMATLAB\CEDS64ML',
                        help='Path to CED MATLAB library')
    parser.add_argument('--cut', action='store_true', default=False,
                        help='Cut strange artifact at the beginning?')
    parser.add_argument('--cut_size', type=float, help='Cut size')

    args = parser.parse_args()

    if args.cut and args.cut_size is None:
        parser.error("--cut requires --cut_size")

    input_dir = args.input_dir
    lfp_freq = args.lfp_freq
    spk_freq = args.spk_freq
    ced_path = args.ceds_lib
    cut = args.cut
    cut_size = args.cut_size

    try:
        load_ceds_lib(ced_path)
    except:
        print 'Cannot load CEDS MATLAB library by path: {}'.format(ced_path)
        USE_MATLAB = False

    walk_res = os.walk(input_dir)

    if USE_MATLAB:
        def write_func(d, r, f): write_with_matlab(d, r, f, spk_freq, lfp_freq)
    else:
        def write_func(d, r, f): write_to_txt(d, r, f)

    for root, _, _ in walk_res:
        map_files = glob.glob(os.path.join(root, '*.map'))
        dat_files = [fname for fname in glob.glob(os.path.join(root, '*.dat')) if os.path.getsize(fname) > 2*BYTES_IN_KB]

        for fname in map_files:
            print 'Working with: {}'.format(fname)

            data_dict = read_alphaomega(fname)
            data_dict = prepare_keys(data_dict)

            write_func(data_dict, root, fname)

        files_by_trial = get_files_by_trials(dat_files)
        for num, files in files_by_trial.items():
            print 'Working with: {}'.format(root)

            data_dict = read_binary(files, cut, cut_size, spk_freq)
            data_dict = prepare_keys(data_dict)

            write_func(data_dict, root, str(num))

    print 'done'
