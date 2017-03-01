import sys

from statistics import *


DEFAULT_BBH = 0.016
DEFAULT_BBL = 0.01
DEFAULT_WIN_SIZE = 13
DEFAULT_REPEAT = 3
DEFAULT_FILT_LOWER = -1
DEFAULT_FILT_UPPER = -1
DEFAULT_FILT_TYPE = FILTER_TIME
DEFAULT_FILENAME = 'burst_res.csv'


PARAMS = dict()

def read_from_file(filename):
    return np.loadtxt(filename)

def get_params():
    global PARAMS
    
    Neuron_Number = 1
    bbh = DEFAULT_BBH
    bbl = DEFAULT_BBL
    brep = DEFAULT_REPEAT
    msize = DEFAULT_WIN_SIZE
    filter_type = DEFAULT_FILT_TYPE
    filter_lower = DEFAULT_FILT_LOWER
    filter_upper = DEFAULT_FILT_UPPER
    file_path = DEFAULT_FILENAME
    data = read_from_file(sys.argv[1])

    if filter_lower == -1:
        if filter_type == FILTER_TIME:
            filter_lower = min(data)
            filter_upper = max(data)
        elif filter_type == FILTER_NUMBER:
            filter_lower = 0
            filter_upper = len(data) + 1



    PARAMS['data'] = data
    PARAMS['bbh'] = bbh
    PARAMS['bbl'] = bbl
    PARAMS['brep'] = brep
    PARAMS['msize'] = msize
    PARAMS['filter'] = filter_type
    PARAMS['filter_lower'] = filter_lower
    PARAMS['filter_upper'] = filter_upper
    PARAMS['file_path'] = file_path + '\\' + DEFAULT_FILENAME

def main():    
    get_params()
    
    spike_data = PARAMS['data']

    time_int = calc_intervals(spike_data, PARAMS['filter'],
                                PARAMS['filter_lower'], PARAMS['filter_upper']) 
                                        
    bi = calc_bi(time_int, PARAMS['bbh'], bbl = PARAMS['bbl'], brep = PARAMS['brep'])
    cv = calc_cv(time_int)
    nu = calc_nu(time_int, PARAMS['msize'])
    br = calc_burst_behavior(time_int)

    print('burst index:', bi)
    print('cv:', cv)
    print('nu:', nu)
    print('burst ratio:', br)
            
            
if __name__ == '__main__':
    main()