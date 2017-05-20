import numpy as np

from nex_info import *

nex_header_dtype_4, nex_variable_dtype_4 = dtype_map[4]

def read_neuron_data_4(f, var_header, file_header):
    f.seek(var_header['data_offset'])
    return np.fromfile(f, dtype=np.int32, count=var_header['count']).astype(np.float64)/file_header['frequency']

def read_event_data_4(f, var_header, file_header):
    f.seek(var_header['data_offset'])
    return np.fromfile(f, dtype=np.int32, count=var_header['count']).astype(np.float64)/file_header['frequency']

def read_interval_data_4(f, var_header, file_header):
    f.seek(var_header['data_offset'])
    begins = np.fromfile(f, dtype=np.int32, count=var_header['count']).astype(np.float64)/file_header['frequency']
    ends = np.fromfile(f, dtype=np.int32, count=var_header['count']).astype(np.float64)/file_header['frequency']
    return np.array(zip(begins, ends))  

read_function_4 = {
    VARIABLE_TYPE_NEURON : read_neuron_data_4,
    VARIABLE_TYPE_EVENT : read_event_data_4,
    VARIABLE_TYPE_INTERVAL : read_interval_data_4,
    # VARIABLE_TYPE_WAVEFORM : read_waveform_data,
    # VARIABLE_TYPE_POPULATION_VECTOR : read_pop_vector_data,
    # VARIABLE_TYPE_CONTINUOUS : read_continuous_data,
    # VARIABLE_TYPE_MARKER : read_marker_data
}
