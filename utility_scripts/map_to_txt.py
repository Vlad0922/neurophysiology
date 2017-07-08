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

OK_SIGNALS_COUNT = 3

def add_quotes(string):
    return '\"{}\"'.format(string)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert AlphaOmega .map files to .txt for spike2')

    parser.add_argument('--input_dir', type=str, required=True,
                        help='Directory with .map files')

    args = parser.parse_args()

    input_dir = args.input_dir

    walk_res = os.walk(input_dir)

    for root, _, _ in walk_res:
        dat_files = glob.glob(os.path.join(root , '*.map'))

        for fname in dat_files:
            print 'Working with: {}'.format(fname)
            r = neo.io.AlphaOmegaIO(filename=fname)

            preproc_path = os.path.join(root, 'converted_to_txt')
            if not os.path.exists(preproc_path):
                os.makedirs(preproc_path)
            
            data_dict = dict()

            for blk in r.read(cascade=True, lazy=False):
                for seg in blk.segments:
                    lfp_signals = sorted(sg.name for sg in seg.analogsignals if sg.name.lower().startswith('lfp'))[:OK_SIGNALS_COUNT]
                    spk_signals = sorted(sg.name for sg in seg.analogsignals if sg.name.lower().startswith('spk'))[:OK_SIGNALS_COUNT]
                    for sg in seg.analogsignals:
                        if sg.name in lfp_signals or sg.name in spk_signals:
                            data_dict[add_quotes(sg.name)] = [v[0] for v in np.array(sg)]

            df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in data_dict.iteritems()]))            
            df = df[map(add_quotes, spk_signals + lfp_signals)]
            df.columns = map(lambda x: x.replace(' ', '_'), df.columns.values)

            fname_no_ext = os.path.splitext(os.path.basename(fname))[0]
            df.to_csv(os.path.join(preproc_path, '{}.txt'.format(fname_no_ext)), index=False, sep='\t')

    print 'done'

