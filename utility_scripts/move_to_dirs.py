# -*- coding: utf-8 -*-

import sys
import os
import glob
import shutil

from collections import defaultdict

BYTES_IN_KB = 1024

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Wrong arguments!')
        print('Usage: python put_in_dirs.py <indir>')
        exit(1)

    input_dir = sys.argv[1]
    walk_res = os.walk(input_dir)
    for root, _, _ in walk_res:
        print 'Working with: {}'.format(root)

        dat_files = [fname for fname in glob.glob(os.path.join(root , '*.dat')) if os.path.getsize(fname) > 2*BYTES_IN_KB]

        files_by_trial = defaultdict(lambda : list())

        for full_path in dat_files:
            fname = os.path.basename(full_path)
            num = fname[:2]
            files_by_trial[num].append(full_path)

        for num, files in files_by_trial.items():
            dir_path = os.path.join(root, num)
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)

            for full_path in files:
                fname = os.path.basename(full_path)
                shutil.copyfile(full_path, os.path.join(dir_path, fname))

    print 'done'