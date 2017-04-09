import numpy as np

from nex_info import *

nex_header_dtype_5, nex_variable_dtype_5 = dtype_map[5]

def read_timestamps_5(f, var_header, file_header):
    if var_header['timestamp_data_type'] == 0:
        return np.fromfile(f, dtype=np.int32, count=var_header['count']).astype(np.float64)/file_header['frequency']
    else:
        return np.fromfile(f, dtype=np.int64, count=var_header['count']).astype(np.float64)/file_header['frequency']


def read_waveform_or_continius_5(f, var_header, file_header):
    points_total = var_header['number_of_data_points']
    if var_header['type'] == VARIABLE_TYPE_WAVEFORM:
        points_total *= var_header['count']

    if var_header['continious_data_type'] == 0:
        return np.fromfile(f, dtype=np.int16, count=points_total).astype(np.float64)*var_header['ad_to_units_coefficient'] + var_header['units_offset']
    else:
        return np.fromfile(f, dtype=np.int32, count=points_total).astype(np.float64)*var_header['ad_to_units_coefficient'] + var_header['units_offset']


def read_neuron_data_5(f, var_header, file_header):
    f.seek(var_header['data_offset'])
    return read_timestamps_5(f, var_header, file_header)

def read_event_data_5(f, var_header, file_header):
    f.seek(var_header['data_offset'])
    return read_timestamps_5(f, var_header, file_header)

def read_interval_data_5(f, var_header, file_header):
    f.seek(var_header['data_offset'])
    starts = read_timestamps_5(f, var_header, file_header)
    ends = read_timestamps_5(f, var_header, file_header)

    return np.array(zip(starts, ends))

def read_waveform_data_5(f, var_header, file_header):
    f.seek(var_header['data_offset'])
    timestamps = read_timestamps_5(f, var_header, file_header)
    mv_data = read_waveform_or_continius_5(f, var_header, file_header)

    res = list()
    for i in range(len(timestamps)):
        wave = list()
        for w in range(var_header['number_of_data_points']):
            wave.append(mv_data[i*var_header['number_of_data_points'] + w])
        res.append(np.array(wave))

    return np.array(wave)

read_function_4 = {
    VARIABLE_TYPE_NEURON : read_neuron_data_5,
    VARIABLE_TYPE_EVENT : read_event_data_5,
    VARIABLE_TYPE_INTERVAL : read_interval_data_5,
    # VARIABLE_TYPE_WAVEFORM : read_waveform_data,
    # VARIABLE_TYPE_POPULATION_VECTOR : read_pop_vector_data,
    # VARIABLE_TYPE_CONTINUOUS : read_continuous_data,
    # VARIABLE_TYPE_MARKER : read_marker_data
}