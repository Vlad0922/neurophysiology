# -*- coding: utf-8 -*-

import sys
import os
import glob

import neo.io

from collections import defaultdict

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import argparse

try:
    import matlab.engine
    eng = matlab.engine.start_matlab()
    eng.addpath('utility_scripts')
    USE_MATLAB = True
except:
    print 'Cannot import or start MATLAB engine, will import to .txt files...'
    USE_MATLAB = False


OK_SIGNALS_COUNT = 3


def add_quotes(string, quotes):
    if quotes:
        return '\"{}\"'.format(string)
    else:
        return string


def load_ceds_lib(lib_path):
    eng.addpath(lib_path, nargout=0)
    eng.CEDS64LoadLib(lib_path, nargout=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert AlphaOmega .map files to .txt/.smr for spike2')

    parser.add_argument('--input_dir', type=str, required=True,
                        help='Directory with .map files')
    parser.add_argument('--with_quotes', action='store_true', default=False,
                        help='Include quotes in header?')
    parser.add_argument('--spk_freq', type=float, required=True,
                        help='Spike channel frequency')
    parser.add_argument('--lfp_freq', type=float, required=True,
                        help='LFP channel frequency')
    parser.add_argument('--ceds_lib', type=str, default=r'C:\CEDMATLAB\CEDS64ML',
                        help='Path to CED MATLAB library')

    args = parser.parse_args()

    input_dir = args.input_dir
    quotes = args.with_quotes
    lfp_freq = args.lfp_freq
    spk_freq = args.spk_freq
    ced_path = args.ceds_lib

    try:
        load_ceds_lib(ced_path)
    except:
        print 'Cannot load CEDS MATLAB library by path: {}'.format(ced_path)
        USE_MATLAB = False

    walk_res = os.walk(input_dir)

    for root, _, _ in walk_res:
        dat_files = glob.glob(os.path.join(root, '*.map'))

        for fname in dat_files:
            print 'Working with: {}'.format(fname)
            r = neo.io.AlphaOmegaIO(filename=fname)

            data_dict = dict()

            for blk in r.read(cascade=True, lazy=False):
                for seg in blk.segments:
                    lfp_signals = sorted(sg.name for sg in seg.analogsignals if sg.name.lower().startswith('lfp'))[:OK_SIGNALS_COUNT]
                    spk_signals = sorted(sg.name for sg in seg.analogsignals if sg.name.lower().startswith('spk'))[:OK_SIGNALS_COUNT]
                    for sg in seg.analogsignals:
                        if sg.name in lfp_signals or sg.name in spk_signals:
                            data_dict[add_quotes(sg.name, quotes)] = [float(v[0]) for v in np.array(sg)]

            fname_no_ext = os.path.splitext(os.path.basename(fname))[0]

            if USE_MATLAB:
                preproc_path = os.path.join(root, 'converted_to_smr')
                if not os.path.exists(preproc_path):
                    os.makedirs(preproc_path)

                for k, v in data_dict.items():
                    data_dict[k] = matlab.double(v)

                orig_keys = data_dict.keys()
                for k in orig_keys:
                    data_dict[k.replace(' ', '_')] = data_dict.pop(k)

                fname_smr = '{}.smr'.format(fname_no_ext)
                fullpath_smr = os.path.join(preproc_path, fname_smr)

                eng.dict_to_smr(data_dict, fullpath_smr, spk_freq, lfp_freq, nargout=0)
            else:
                preproc_path = os.path.join(root, 'converted_to_txt')
                if not os.path.exists(preproc_path):
                    os.makedirs(preproc_path)

                df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in data_dict.iteritems()]))
                df = df[map(lambda x: add_quotes(x, quotes), spk_signals + lfp_signals)]
                df.columns = map(lambda x: x.replace(' ', '_'), df.columns.values)

                df.to_csv(os.path.join(preproc_path, '{}.txt'.format(fname_no_ext)), index=False, sep='\t')

            exit(1)

    print 'done'
