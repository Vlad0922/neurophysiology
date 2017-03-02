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
    file_path = os.environ['HOMEPATH'] + '\\Desktop\\'
    freq = DEFAULT_FREQ
    fOscillationFreq = 35
    Fmin = 30
    Fmax = 50
    interval = 1
    
    doc = nex.GetActiveDocument()
    __wrapper = []
    res = nex.Dialog(doc, 
                    Neuron_Number, "Select Neuron", "neuron", 
                    interval, "Select interval(s)", "interval",
                    bbh, "Upper bound", "number",
                    bbl, "Lower bound", "number",
                    brep, "Repetition count", "number",
                    msize, "Window size", "number",
					freq, "Signal frequency", "number",
					fOscillationFreq, "Oscillation frequency", "number",
					Fmin, "Min filter frequency", "number",
					Fmax, "Max filter frequency", "number",
                    file_path, "Result file path", "string", 
                    __wrapper )       
    
    Neuron_Number = __wrapper[0]
    Neuron_Var = nex.GetVar(doc, Neuron_Number, "neuron")
    
    interval = __wrapper[1]
    interval_var = nex.GetVar(doc, interval, "interval")
    
    bbh = float(__wrapper[2])
    bbl = float(__wrapper[3])
    brep = int(__wrapper[4])
    msize = int(__wrapper[5])
    freq = int(__wrapper[6])
    fOscillationFreq = int(__wrapper[7])
    Fmin = int(__wrapper[8])
    Fmax = int(__wrapper[9])
    file_path = str(__wrapper[10])
    
    PARAMS['data'] = np.array(Neuron_Var.Timestamps())
    PARAMS['bbh'] = bbh
    PARAMS['bbl'] = bbl
    PARAMS['brep'] = brep
    PARAMS['msize'] = msize
    PARAMS['filters'] = [p for p in zip(interval_var.Intervals()[0], interval_var.Intervals()[1])]
    PARAMS['data_name'] = nex.GetName(Neuron_Var)
    PARAMS['file_path'] = file_path + '\\' + DEFAULT_FILENAME
    PARAMS['frequency'] = freq
    PARAMS['Fmin'] = Fmin
    PARAMS['Fmax'] = Fmax


def main():    
    get_params()
    
    data_raw = PARAMS['data']
    data_filtered = np.array([])
    time_int = np.array([])
    for filter in PARAMS['filters']:
        spike_data = filter_spikes(data_raw, filter[0], filter[1])
        data_filtered = np.concatenate([data_filtered, spike_data])
        time_int = np.concatenate([time_int, calc_intervals(spike_data)])
            
    if len(time_int) == 0:
        raise 'Empty filter result!'
                                        
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

    Trial = sec_to_timestamps(data_filtered, PARAMS['frequency']).tolist()
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