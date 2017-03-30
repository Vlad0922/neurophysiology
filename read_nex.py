import numpy as np

from nex_info import *

def read_timestamps(f, var_header, file_header):
    if var_header['timestamp_data_type'] == 0:
        return np.fromfile(f, dtype=np.int32, count=var_header['count']).astype(np.float64)/file_header['frequency']
    else:
        return np.fromfile(f, dtype=np.int64, count=var_header['count']).astype(np.float64)/file_header['frequency']


def read_waveform_or_continius(f, var_header, file_header):
    points_total = var_header['number_of_data_points']
    if var_header['type'] == VARIABLE_TYPE_WAVEFORM:
        points_total *= var_header['count']

    if var_header['continious_data_type'] == 0:
        return np.fromfile(f, dtype=np.int16, count=points_total).astype(np.float64)*var_header['ad_to_units_coefficient'] + var_header['units_offset']
    else:
        return np.fromfile(f, dtype=np.int32, count=points_total).astype(np.float64)*var_header['ad_to_units_coefficient'] + var_header['units_offset']


def read_neuron_data(f, var_header, file_header):
    f.seek(var_header['data_offset'])
    return read_timestamps(f, var_header, file_header)

def read_event_data(f, var_header, file_header):
    f.seek(var_header['data_offset'])
    return read_timestamps(f, var_header, file_header)

def read_interval_data(f, var_header, file_header):
    f.seek(var_header['data_offset'])
    starts = read_timestamps(f, var_header, file_header)
    ends = read_timestamps(f, var_header, file_header)

    return np.array(zip(starts, ends))

def read_waveform_data(f, var_header, file_header):
    f.seek(var_header['data_offset'])
    timestamps = read_timestamps(f, var_header, file_header)
    mv_data = read_waveform_or_continius(f, var_header, file_header)

    res = list()
    for i in range(len(timestamps)):
        wave = list()
        for w in range(var_header['number_of_data_points']):
            wave.append(mv_data[i*var_header['number_of_data_points'] + w])
        res.append(np.array(wave))

    return np.array(wave)


read_function = {
    VARIABLE_TYPE_NEURON : read_neuron_data,
    VARIABLE_TYPE_EVENT : read_event_data,
    VARIABLE_TYPE_INTERVAL : read_interval_data,
    # VARIABLE_TYPE_WAVEFORM : read_waveform_data,
    # VARIABLE_TYPE_POPULATION_VECTOR : read_pop_vector_data,
    # VARIABLE_TYPE_CONTINUOUS : read_continuous_data,
    # VARIABLE_TYPE_MARKER : read_marker_data
}

with open('TestDataFileForHowTo.nex5', 'r') as f:
    nex_header = np.fromfile(f, dtype=nex_header_dtype, count=1)

    vars_info = np.fromfile(f, dtype=nex_variable_dtype, count=nex_header['number_of_variables'])
    vars_data = list()

    for v in vars_info:
        var_type = v['type']
        if var_type in read_function:
            vars_data.append(read_function[var_type](f, v, nex_header))
            print v['name'], vars_data[-1][:5], len(vars_data[-1])
