import numpy as np

from nex_5_reader import *
from nex_4_reader import *

def get_version(f):
    magic = np.fromfile(f, dtype=np.int32, count=1)[0]
    f.seek(0)
    if magic == NEX4_MAGIC:
        return 4
    elif magic == NEX5_MAGIC:
        return 5
    else:
        return -1

def read_nex_file(fname):
    with open(fname, 'rb') as f:
        nex_ver = get_version(f)

        print 'has nex{} file format'.format(nex_ver)
        if nex_ver == 5:
            read_function = read_function_5
            nex_header_dtype, nex_variable_dtype = nex_header_dtype_5, nex_variable_dtype_5
        elif nex_ver == 4:
            read_function = read_function_4
            nex_header_dtype, nex_variable_dtype = nex_header_dtype_4, nex_variable_dtype_4
        
        nex_header = np.fromfile(f, dtype=nex_header_dtype, count=1)

        vars_info = np.fromfile(f, dtype=nex_variable_dtype, count=nex_header['number_of_variables'])
        vars_data = list()

        for v in vars_info:
            var_type = v['type']
            if var_type in read_function:
                vars_data.append((v['name'], var_type, read_function[var_type](f, v, nex_header)))

        return vars_data
                
