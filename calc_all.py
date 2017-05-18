from __future__ import division

import os

import utility
reload(utility)
from utility import *

import statistics
reload(statistics)
from statistics import *

import oscore
reload(oscore)
from oscore import *

import nex

DEFAULT_FILENAME = 'result.xlsx'
PARAMS = dict()   
ISP_RANGE = [(1,3), (3,8), (8,13), (13,30), (30,100)]
#OSCORE_RANGE = [(1., 3.), (3.,8.), (8.,12.), (12.,20.), (20.,30.), (30.,60.), (60., 90.)]
OSCORE_RANGE = [(3.,8.), (8.,12.), (12.,20.), (20.,30.), (30.,60.), (60., 90.)]

def get_params():
    global PARAMS
    
    Neuron_Number = 1
    bbh = DEFAULT_BBH
    bbl = DEFAULT_BBL
    brep = DEFAULT_REPEAT
    msize = DEFAULT_WIN_SIZE
    file_path = os.environ['HOMEPATH'] + '\\Desktop\\'
    freq = DEFAULT_FREQ
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
                    freq, "Ms per bin", "number",
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
    freq = float(__wrapper[6])
    Fmin = int(__wrapper[7])
    Fmax = int(__wrapper[8])
    file_path = str(__wrapper[9])
    
    PARAMS['data'] = np.array(Neuron_Var.Timestamps())
    PARAMS['bbh'] = bbh
    PARAMS['bbl'] = bbl
    PARAMS['brep'] = brep
    PARAMS['msize'] = msize
    PARAMS['filters'] = [p for p in zip(interval_var.Intervals()[0], interval_var.Intervals()[1])]
    PARAMS['data_name'] = nex.GetName(Neuron_Var)
    PARAMS['file_path'] = file_path + '\\' + DEFAULT_FILENAME
    PARAMS['frequency'] = int(1/freq)
    PARAMS['Fmin'] = Fmin
    PARAMS['Fmax'] = Fmax
    PARAMS['doc_name'] = nex.GetDocTitle(doc)


def main():    
    get_params()
    
    data_raw = PARAMS['data']
    data_filtered = np.array([])
    time_int = np.array([])
    interval_len = 0.
    for l, r in PARAMS['filters']:
        interval_len += r - l
         
        spike_data = filter_spikes(data_raw, l, r)
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

#    Trial = sec_to_timestamps(data_filtered, PARAMS['frequency']).tolist()
#    iTrialLength = Trial[-1]
#    oscore = oscore_spikes(np.array([Trial]), iTrialLength, PARAMS['Fmin'], PARAMS['Fmax'], PARAMS['frequency'])
    
    df = dict()
    df['data_name'] = PARAMS['data_name']
    df['burst index'] = bi
    df['cv'] = cv
    df['nu'] = nu
    df['frequency variance'] = freq_v
    df['modalirity burst'] = mod_burst
    df['pause index'] = pause_ind
    df['pause ratio'] = pause_rat
    df['skewness'] = skew
    df['kurtoisis'] = kurt
    df['burst_mean'] = burst_mean
#    df['oscore'] = oscore
    df['type'] = get_type(df['burst_mean'], df['cv'])
    df['doc_name'] = PARAMS['doc_name']
    df['isi_mean'] = np.mean(time_int)
    df['isi_median'] = np.median(time_int)
    df['isi_std'] = np.std(time_int)
    df['spike_count'] = len(data_filtered)
    df['filter_length'] = interval_len
    df['bi_2'] = calc_bi_two(data_filtered)
    df['lv'] = calc_local_variance(time_int)
    df['firing_rate'] = 1.*(data_filtered[-1] - data_filtered[0])/len(data_filtered)
    df['burst_rate'] = calc_burst_rate(time_int)
    
    for (osc_l, osc_h) in OSCORE_RANGE:
        Trial = sec_to_timestamps(data_filtered, PARAMS['frequency']).tolist()
        iTrialLength = Trial[-1]
        oscore = oscore_spikes(np.array([Trial]), iTrialLength, osc_l, osc_h, PARAMS['frequency'])
        df['oscore_{}_{}'.format(osc_l, osc_h)] = oscore

    #for (isp_l, isp_h) in ISP_RANGE:
    #    df['ISP_{}_{}'.format(isp_l, isp_h)] = calc_isp(time_int, isp_l, isp_h)
    
    write_to_excel(PARAMS['file_path'], 'all_results', df, ['doc_name', 'data_name'])
    for key in df.keys():
        df[key] = [df[key]]
    
    #draw_table(df)
    print 'done'
            
            
if __name__ == '__main__':
    main()