import sys
import os

from read_nex import read_nex_file

from statistics import *
from utility import *


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print 'Wrong usage!\n'
        print 'use python handle_files.py input_dir output_file\n'
        exit(1)

    input_path = sys.argv[1]
    output_fname = sys.argv[2]

    for file in os.listdir(input_path):
        file_path = input_path + file
        sys.stdout.write('reading file {}...'.format(file_path))
        file_data = read_nex_file(file_path)
        print len(file_data)
