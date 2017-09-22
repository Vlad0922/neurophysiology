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
                    wrapper )   
    
    phase_var = nex.GetVar(doc, wrapper[0], "continuous")
    amp_var = nex.GetVar(doc, wrapper[1], "continuous")
    
    params['phase_data'] = np.array(phase_var.Timestamps())
    params['amp_data'] = np.array(amp_var.Timestamps())
    params['f_range_lo'] = map(float, wrapper[2].split(','))
    params['f_range_hi'] = map(float, str(wrapper[3]).split(','))
    params['seconds_lo'] = float(wrapper[4])
    params['seconds_hi'] = float(wrapper[5])
    params['fs'] = phase_var.SamplingRate()
    params['method'] = str(wrapper[6])
    
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
    print pac


main()