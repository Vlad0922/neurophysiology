import numpy as np


VARIABLE_TYPE_NEURON    = 0
VARIABLE_TYPE_EVENT     = 1
VARIABLE_TYPE_INTERVAL  = 2
VARIABLE_TYPE_WAVEFORM  = 3
VARIABLE_TYPE_POPULATION_VECTOR = 4
VARIABLE_TYPE_CONTINUOUS        = 5
VARIABLE_TYPE_MARKER            = 6

nex_header_dtype = np.dtype([
    ('magic_number', '<i4'),
    ('file_version', '<i4'),
    ('comment', 'S256'),
    ('frequency', np.float64),
    ('start_time_in_ticks', '<i8'),
    ('number_of_variables', '<i4'),
    ('metadata_offset', '<u8'),
    ('padding', (bytes, 64))
    ], align=0)


nex_variable_dtype = np.dtype([
    ('type', np.int32),
    ('version', np.int32),
    ('name', (bytes, 64)),
    ('data_offset', np.uint64),
    ('count', np.uint64),
    ('timestamp_data_type', '<i4'),
    ('continious_data_type', '<i4'),
    ('frequency', np.float64),
    ('units', (bytes, 32)),
    ('ad_to_units_coefficient', np.float64),
    ('units_offset', np.float64),
    ('number_of_data_points', '<u8'),
    ('prethreshold_time', np.float64),
    ('marker_data_type', '<i4'),
    ('number_of_marker_fields', '<i4'),
    ('marker_len', '<i4'),
    ('index_of_first_point', '<i4'),
    ('padding', (bytes, 60))
    ], align=0)