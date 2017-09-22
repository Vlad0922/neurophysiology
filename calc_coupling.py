import nex
import sys

import matplotlib.pyplot as plt

import numpy as np
import neurodsp
#from pacpy.pac import plv

def get_params():
    params = dict()
        
    phase_neuron = 1
    amp_neuron = 2
    f_range_lo = '13,30'
    f_range_hi = '50,200'
    seconds_lo = 0.25
    seconds_hi = 0.2
    method = 'ozkurt'
    secs = 50
    f_pha_bin_edges = '2,42,2' #np.arange(2, 42, 2)
    f_amp_bin_edges = '2,200,4' #np.arange(20, 200, 4)
    N_cycles_pha = 5
    N_cycles_amp = 11
    

    doc = nex.GetActiveDocument()
    wrapper = []
    res = nex.Dialog(doc, 
                    phase_neuron, "Neuron from which to compute the phase component", "continuous", 
                    amp_neuron, "Neuron from which to compute the amplitude component", "continuous",
                    f_range_lo, "The low frequency filtering range (Hz)", "string",
                    f_range_hi, "The high frequency filtering range (Hz)", "string",
                    seconds_lo, "Length of the low band-pass filter (seconds)", "number",
                    seconds_hi, "Length of the high band-pass filter (seconds)", "number",
                    method, "Method to compute the PAC", "string",
                    secs, "Length of the interval to compute the PAC", "number",
                    f_pha_bin_edges, 'Min, max and step for array of edges of the freq.r. for phase', 'string',
                    f_amp_bin_edges, 'Min, max and step for array of edges of the freq.r. for amplitude', 'string',
                    N_cycles_pha, 'Cycles count for the phase filter', 'number',
                    N_cycles_amp, 'Cycles count for the amplitude filter', 'number',
                    wrapper )   
    
    phase_var = nex.GetVar(doc, wrapper[0], "continuous")
    amp_var = nex.GetVar(doc, wrapper[1], "continuous")
    
    secs = int(wrapper[7])
    idx = np.where(np.array(phase_var.Timestamps()) <= secs)
        
    params['fs'] = phase_var.SamplingRate()    
    params['phase_data'] = np.array(phase_var.ContinuousValues())[idx]
    params['phase_name'] = phase_var.Name()
    params['amp_data'] = np.array(amp_var.ContinuousValues())[idx]
    params['amp_name'] = amp_var.Name()
    params['f_range_lo'] = map(float, wrapper[2].split(','))
    params['f_range_hi'] = map(float, str(wrapper[3]).split(','))
    params['seconds_lo'] = float(wrapper[4])
    params['seconds_hi'] = float(wrapper[5])    
    params['method'] = str(wrapper[6])
    
    pha_b, pha_e, pha_s = map(float, str(wrapper[8]).split(','))
    amp_b, amp_e, amp_s = map(float, str(wrapper[9]).split(','))
    
    params['f_pha_bin_edges'] = np.arange(pha_b, pha_e, pha_s)
    params['f_amp_bin_edges'] = np.arange(amp_b, amp_e, amp_s)
    params['N_cycles_pha'] = int(wrapper[10])
    params['N_cycles_amp'] = int(wrapper[11])
    
    return params


def main():
    params = get_params()
    pac = neurodsp.compute_pac(params['phase_data'],
                               params['amp_data'],
                               params['fs'], 
                               params['f_range_lo'],
                               params['f_range_hi'],
                               N_seconds_lo=params['seconds_lo'],
                               N_seconds_hi=params['seconds_hi'],
                               pac_method=params['method'])
    print 'PAC={}'.format(pac)
    
#    t = np.arange(0, params['phase_data'].shape[0]/params['fs'], 1/params['fs'])
#    plt.plot(t, params['phase_data'])
#    plt.show()
    
    pac_plot = neurodsp.compute_pac_comodulogram(params['phase_data'],
                                                    params['amp_data'],
                                                    params['fs'],
                                                    params['f_pha_bin_edges'],
                                                    params['f_amp_bin_edges'],
                                                    N_cycles_pha=params['N_cycles_pha'],
                                                    N_cycles_amp=params['N_cycles_amp'],
                                                    pac_method=params['method'])
                                        
    neurodsp.plot_pac_comodulogram(pac_plot, params['f_pha_bin_edges'], params['f_amp_bin_edges'],
                               clim=(0,.1), colormap='Spectral')
    ax = plt.gca()
    xlabel = ax.get_xlabel()
    ylabel = ax.get_ylabel() 
    
    ax.set_xlabel(xlabel + ' ' + params['phase_name'])
    ax.set_ylabel(ylabel + ' ' + params['amp_name'])
    
    plt.title('PAC value={}'.format(pac))
    plt.show()
    

main()