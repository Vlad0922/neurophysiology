# -*- coding: utf-8 -*-

import sys
import os
import glob
import shutil

from collections import defaultdict

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import argparse

BYTES_IN_KB = 1024
UV_IN_VOLTS = 0.001
FREQ = 20000

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert binary files to .txt for spike2')

    parser.add_argument('--input_dir', type=str, required=True,
                        help='Directory with binary files')
    parser.add_argument('--cut', action='store_true', default=True,
                        help='Cut strange artifact at the beginning?')

    args = parser.parse_args()

    input_dir = args.input_dir
    cut = args.cut
    walk_res = os.walk(input_dir)

    for root, _, _ in walk_res:
        dat_files = [fname for fname in glob.glob(os.path.join(root , '*.dat')) if os.path.getsize(fname) > 2*BYTES_IN_KB]

        files_by_trial = defaultdict(lambda : list())

        for full_path in dat_files:
            fname = os.path.basename(full_path)
            num = fname[:2]
            files_by_trial[num].append(full_path)

        if len(files_by_trial) != 0:
            print 'Working with: {}'.format(root)
            preproc_path = os.path.join(root, 'converted_to_txt')
            if not os.path.exists(preproc_path):
                os.makedirs(preproc_path)

            for num, files in files_by_trial.items():
                data_dict = dict()
                for full_path in files:
                    fname = os.path.basename(full_path)
                    data = np.fromfile(full_path, dtype="<u2")*UV_IN_VOLTS

                    if cut:
                        data = data[int(0.2*FREQ):]

                    data_dict['\"{}\"'.format(fname)] = data

                orig_cols = data_dict.keys()

                max_len = np.max([len(data) for data in data_dict.values()])
                data_dict['"Time"'] = np.arange(max_len, dtype=float)/FREQ

                df = pd.DataFrame(data=data_dict)
                df = df[['"Time"'] + orig_cols]

                df.to_csv(os.path.join(preproc_path, '{}.txt'.format(str(num))), index=False, sep='\t')
    print 'done'