# -*- coding: utf-8 -*-

import sys
import os
import glob

from collections import defaultdict

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import argparse


BYTES_IN_KB = 1024
UV_IN_VOLTS = 0.001
FREQ = 20000.
DUMMY = FREQ


try:
    import matlab.engine
    eng = matlab.engine.start_matlab()
    eng.addpath('utility_scripts')
    USE_MATLAB = True
except:
    print 'Cannot import or start MATLAB engine, will import to .txt files...'
    USE_MATLAB = False


def load_ceds_lib(lib_path):
    eng.addpath(lib_path, nargout=0)
    eng.CEDS64LoadLib(lib_path, nargout=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert binary files to .txt for spike2')

    parser.add_argument('--input_dir', type=str, required=True,
                        help='Directory with binary files')
    parser.add_argument('--cut', action='store_true', default=False,
                        help='Cut strange artifact at the beginning?')
    parser.add_argument('--ceds_lib', type=str, default=r'C:\CEDMATLAB\CEDS64ML',
                    help='Path to CED MATLAB library')
    parser.add_argument('--with_matlab', action='store_true', default=False,
                        help='Use MATLAB for convertation?')

    args = parser.parse_args()

    input_dir = args.input_dir
    cut = args.cut
    ced_path = args.ceds_lib

    try:
        load_ceds_lib(ced_path)
    except:
        print 'Cannot load CEDS MATLAB library by path: {}'.format(ced_path)
        USE_MATLAB = False

    USE_MATLAB &= args.with_matlab

    walk_res = os.walk(input_dir)

    try:
        load_ceds_lib(args.ceds_lib)
    except:
        print 'Cannot load CEDS MATLAB library by path: {}'.format(args.ceds_lib)
        USE_MATLAB = False

    for root, _, _ in walk_res:
        dat_files = [fname for fname in glob.glob(os.path.join(root , '*.dat')) if os.path.getsize(fname) > 2*BYTES_IN_KB]

        files_by_trial = defaultdict(lambda : list())

        for full_path in dat_files:
            fname = os.path.basename(full_path)
            num = fname[:2]
            files_by_trial[num].append(full_path)

        for num, files in files_by_trial.items():
            print 'Working with: {}'.format(root)

            data_dict = dict()
            for full_path in files:
                fname = os.path.splitext(os.path.basename(full_path))[0]
                fname = fname[fname.find('_')+1:]

                data = np.fromfile(full_path, dtype="<u2")*UV_IN_VOLTS

                if cut:
                    data = data[int(0.2*FREQ):]

                data_dict[fname] = data

            if USE_MATLAB:
                preproc_path = os.path.join(root, 'converted_to_smr')
                if not os.path.exists(preproc_path):
                    os.makedirs(preproc_path)

                for k, v in data_dict.items():
                    data_dict[k] = matlab.double([float(val) for val in v])

                orig_keys = data_dict.keys()
                for k in orig_keys:
                    data_dict[k.replace(' ', '_')] = data_dict.pop(k)

                fname_smr = '{}.smr'.format(str(num))
                fullpath_smr = os.path.join(preproc_path, fname_smr)

                eng.dict_to_smr(data_dict, fullpath_smr, FREQ, DUMMY, nargout=0)
            else:
                preproc_path = os.path.join(root, 'converted_to_txt')
                if not os.path.exists(preproc_path):
                    os.makedirs(preproc_path)
                orig_cols = data_dict.keys()

                max_len = np.max([len(data) for data in data_dict.values()])
                data_dict['Time'] = np.arange(max_len, dtype=float)/FREQ

                df = pd.DataFrame(data=data_dict)
                df = df[['Time'] + orig_cols]

                df.to_csv(os.path.join(preproc_path, '{}.txt'.format(str(num))), index=False, sep='\t')

            exit(1)

    print 'done'