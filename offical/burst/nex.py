import os.path

import nex

DEFAULT_FILENAME = 'burst_res.csv'


PARAMS = dict()


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
    
    doc = nex.GetActiveDocument()
    __wrapper = []
    res = nex.Dialog(doc, 
                    Neuron_Number, "Select Neuron", "neuron", 
                    bbh, "Upper bound", "number",
                    bbl, "Lower bound", "number",
                    brep, "Repetition count", "number",
                    msize, "Window size", "number",
                    filter_type, "Filter by time or interval number", "number",
                    filter_lower, "Lower bound of filtering", "number",
                    filter_upper, "Upper bound of filtering", "number",
                    file_path, "Result file path", "string", 
                    __wrapper )       
    
    Neuron_Number = __wrapper[0];
    Neuron_Var = nex.GetVar(doc, Neuron_Number, "neuron")
    
    bbh = float(__wrapper[1])
    bbl = float(__wrapper[2])
    brep = int(__wrapper[3])
    msize = int(__wrapper[4])
    filter_type = int(__wrapper[5])
    filter_lower = float(__wrapper[6])
    filter_upper = float(__wrapper[7])
    file_path = str(__wrapper[8])
    
    if filter_lower == -1:
        if filter_type == 1:
            filter_lower = min(Neuron_Var.Timestamps())
            filter_upper = max(Neuron_Var.Timestamps())
        elif filter_type == 0:
            filter_lower = 0
            filter_upper = len(Neuron_Var.Timestamps()) + 1
    
    PARAMS['data'] = np.array(Neuron_Var.Timestamps())
    PARAMS['bbh'] = bbh
    PARAMS['bbl'] = bbl
    PARAMS['brep'] = brep
    PARAMS['msize'] = msize
    PARAMS['filter'] = filter_type
    PARAMS['filter_lower'] = filter_lower
    PARAMS['filter_upper'] = filter_upper
    PARAMS['data_name'] = nex.GetName(Neuron_Var)
    PARAMS['file_path'] = file_path + '\\' + DEFAULT_FILENAME

def main():    
    get_params()
    
    spike_data = PARAMS['data']
        
    time_int = calc_intervals(spike_data, PARAMS['filter'],
                                PARAMS['filter_lower'], PARAMS['filter_upper']) 
                                        
    bi = calc_bi(time_int, PARAMS['bbh'], bbl = PARAMS['bbl'], brep = PARAMS['brep'])
    cv, nu = calc_cv(time_int)
    
    if os.path.isfile(PARAMS['file_path']):
        with open(PARAMS['file_path'], 'a+') as out:
            out.write('{};{};{};{}\n'.format(PARAMS['data_name'], bi, cv, nu))
    else:
        with open(PARAMS['file_path'], 'w') as out:
            out.write('data_name;bi;cv;nu\n')
            out.write('{};{};{};{}\n'.format(PARAMS['data_name'], bi, cv, nu))
            
    print 'done'
            
            
if __name__ == '__main__':
    main()