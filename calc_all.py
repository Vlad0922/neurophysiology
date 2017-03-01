from __future__ import division

import os.path
import imp

from statistics import *

from oscore import *

import nex

DEFAULT_FILENAME = 'all_res.csv'

PARAMS = dict()


def sec_to_timestamps(spikes, freq):
    return np.array((spikes-np.floor(spikes[0]))*PARAMS['frequency'], dtype=np.int32)
	

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
    file_path = os.environ['HOMEPATH'] + '\\Desktop\\'
    freq = DEFAULT_FREQ
    fOscillationFreq = 35
    Fmin = 30
    Fmax = 50
    
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
					freq, "Signal frequency", "number",
					fOscillationFreq, "Oscillation frequency", "number",
					Fmin, "Min filter frequency", "number",
					Fmax, "Max filter frequency", "number",
                    file_path, "Result file path", "string", 
                    __wrapper )       
    
    Neuron_Number = __wrapper[0]
    Neuron_Var = nex.GetVar(doc, Neuron_Number, "neuron")
    
    bbh = float(__wrapper[1])
    bbl = float(__wrapper[2])
    brep = int(__wrapper[3])
    msize = int(__wrapper[4])
    filter_type = int(__wrapper[5])
    filter_lower = float(__wrapper[6])
    filter_upper = float(__wrapper[7])
    freq = int(__wrapper[8])
    fOscillationFreq = int(__wrapper[9])
    Fmin = int(__wrapper[10])
    Fmax = int(__wrapper[11])
    file_path = str(__wrapper[12])
    
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
    PARAMS['frequency'] = freq
    PARAMS['Fmin'] = Fmin
    PARAMS['Fmax'] = Fmax


def main():    
    get_params()
    
    spike_data = PARAMS['data']
        
    time_int = calc_intervals(spike_data, filter_by=PARAMS['filter'],
                                low=PARAMS['filter_lower'], high=PARAMS['filter_upper']) 
                                        
    bi = calc_burst_index(time_int, bbh=PARAMS['bbh'], bbl=PARAMS['bbl'], brep=PARAMS['brep'])
    cv = calc_cv(time_int)
    nu = calc_nu(time_int)
    freq_v = calc_freq_var(time_int)
    mod_burst = calc_modalirity_burst(time_int)
    pause_ind = calc_pause_index(time_int)
    pause_rat = calc_pause_ratio(time_int)
    burst_beh = calc_burst_behavior(time_int)
    skew = calc_skewness(time_int)
    kurt = calc_kurtosis(time_int)
    burst_mean = calc_burst_by_mean(time_int)

    Trial = sec_to_timestamps(PARAMS['data'], PARAMS['frequency']).tolist()
    iTrialLength = Trial[-1]
    oscore = oscore_spikes(np.array([Trial]), iTrialLength, PARAMS['Fmin'], PARAMS['Fmax'], PARAMS['frequency'])
    
    if not os.path.isfile(PARAMS['file_path']):
        with open(PARAMS['file_path'], 'w') as out:
            out.write('data_name,burst_index,cv,nu,frequence_variance,modalirity_burst,'
                        'pause_index,pause_ratio,burst_behavior,skewness,kurtosis,oscore,burst_mean\n')
        
    with open(PARAMS['file_path'], 'a+') as out:
            out.write('{},{},{},{},{},{},{},{},{},{},{},{},{}\n'.format(PARAMS['data_name'], bi, 
                                             cv, nu, freq_v, mod_burst, pause_ind,
                                             pause_rat, burst_beh, skew, kurt, oscore, burst_mean))
            
    print 'done'
            
            
if __name__ == '__main__':
    main()