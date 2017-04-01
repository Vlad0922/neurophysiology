import numpy as np

NEX4_MAGIC = 827868494
NEX5_MAGIC = 894977358

VARIABLE_TYPE_NEURON    = 0
VARIABLE_TYPE_EVENT     = 1
VARIABLE_TYPE_INTERVAL  = 2
VARIABLE_TYPE_WAVEFORM  = 3
VARIABLE_TYPE_POPULATION_VECTOR = 4
VARIABLE_TYPE_CONTINUOUS        = 5
VARIABLE_TYPE_MARKER            = 6

nex5_header_dtype = np.dtype([
    ('magic_number', np.int32),
    ('file_version', np.int32),
    ('comment', 'S256'),
    ('frequency', np.float64),
    ('start_time_in_ticks', np.int64),
    ('number_of_variables', np.int32),
    ('metadata_offset', np.uint64),
    ('padding', (bytes, 64))
    ], align=0)


nex5_variable_dtype = np.dtype([
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


nex4_header_dtype = np.dtype([
    ('magic_number', np.int32),
    ('file_version', np.int32),
    ('comment', 'S256'),
    ('frequency', np.float64),
    ('begin', np.int32),
    ('end', np.int32),
    ('number_of_variables', np.int32),
    ('next_file_header', np.int32),
    ('padding', (bytes, 256))
    ], align=0)


nex4_variable_dtype = np.dtype([
    ('type', np.int32),
    ('version', np.int32),
    ('name', 'S64'),
    ('data_offset', np.int32),
    ('count', np.int32),
    ('wire_num', np.int32),
    ('unit_num', np.int32),
    ('gain', np.int32),
    ('filter', np.int32),
    ('xpos', np.float64),
    ('ypos', np.float64),
    ('wfrequency', np.float64),
    ('adtomv', np.float64),
    ('npoints_wave', np.int32),
    ('nmarkers', np.int32),
    ('marker_len', np.int32),
    ('mv_offset', np.float64),
    ('prethreshold_time', np.float64),
    ('padding', (bytes, 52))
    ])


dtype_map = {
    5 : (nex5_header_dtype, nex5_variable_dtype),
    4 : (nex4_header_dtype, nex4_variable_dtype)
}