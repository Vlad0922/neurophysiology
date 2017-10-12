import nex
import sys

import matplotlib.pyplot as plt

import numpy as np
import neurodsp


def get_params():
    params = dict()
        
    phase_neuron = 1
    secs = 50.

    doc = nex.GetActiveDocument()
    wrapper = []
    res = nex.Dialog(doc, 
                    phase_neuron, "Neuron from which to compute the phase component", "continuous", 
                    secs,  "Length of the interval to compute the PSD", "number",
                    wrapper )   
    
    phase_var = nex.GetVar(doc, wrapper[0], "continuous")
    
    secs = int(wrapper[1])
    idx = np.where(np.array(phase_var.Timestamps()) <= secs)
        
    params['fs'] = phase_var.SamplingRate()    
    params['phase_data'] = np.array(phase_var.ContinuousValues())[idx]
    params['var_name'] = phase_var.Name()
    
    return params
    

def main():
    params = get_params()
    
    x = params['phase_data']
    fs = params['fs']

    freq_mean, P_mean = neurodsp.spectral.psd(x, fs, method='mean', nperseg=fs*2) # mean of spectrogram (Welch)
    freq_med, P_med = neurodsp.spectral.psd(x, fs, method='median', nperseg=fs*2) # median of spectrogram ("median Welch")
    freq_mf, P_mf = neurodsp.spectral.psd(x, fs, method='medfilt')

    plt.figure(figsize=(8,8))
    plt.loglog(freq_mean[:200],P_mean[:200], label='Welch')
    plt.loglog(freq_med[:200],P_med[:200], label='Median Welch')
    plt.loglog(freq_mf[100:10000],P_mf[100:10000], label='Median Filter FFT')
    plt.legend()
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power (V^2/Hz)')
    plt.title('PSD of {}'.format(params['var_name']))    
    plt.show()   
    
    slope, offset= neurodsp.spectral.fit_slope(freq_med[1:200], P_med[1:200], (30., 100.), fit_excl=[(55., 65.)], method='ols', plot_fit=True)
    plt.title('1/f noise of {}'.format(params['var_name']))
    plt.show()
    
    
main()