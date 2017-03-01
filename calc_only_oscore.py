from __future__ import division

import os.path
import imp

from oscore import *

import nex

DEFAULT_FILENAME = 'oscore_res.csv'

PARAMS = dict()


def sec_to_timestamps(spikes, freq):
    return np.array((spikes-np.floor(spikes[0]))*PARAMS['frequency'], dtype=np.int32)
	

def get_params():
    global PARAMS
    
    Neuron_Number = 1
    file_path = os.environ['HOMEPATH'] + '\\Desktop\\'
    freq = DEFAULT_FREQ
    fOscillationFreq = 35
    Fmin = 30
    Fmax = 50
    
    doc = nex.GetActiveDocument()
    __wrapper = []
    res = nex.Dialog(doc, 
                    Neuron_Number, "Select Neuron", "neuron", 
					freq, "Signal frequency", "number",
					fOscillationFreq, "Oscillation frequency", "number",
					Fmin, "Min filter frequency", "number",
					Fmax, "Max filter frequency", "number",
                    file_path, "Result file path", "string", 
                    __wrapper )       
    
    Neuron_Number = __wrapper[0]
    Neuron_Var = nex.GetVar(doc, Neuron_Number, "neuron")
    
    freq = int(__wrapper[1])
    fOscillationFreq = int(__wrapper[2])
    Fmin = int(__wrapper[3])
    Fmax = int(__wrapper[4])
    file_path = str(__wrapper[5])
    
    
    PARAMS['data'] = np.array(Neuron_Var.Timestamps())
    PARAMS['data_name'] = nex.GetName(Neuron_Var)
    PARAMS['file_path'] = file_path + '\\' + DEFAULT_FILENAME
    PARAMS['frequency'] = freq
    PARAMS['Fmin'] = Fmin
    PARAMS['Fmax'] = Fmax


def main():    
    get_params()
    
    spike_data = PARAMS['data']

    Trial = sec_to_timestamps(PARAMS['data'], PARAMS['frequency']).tolist()
    iTrialLength = Trial[-1]
    oscore = oscore_spikes(np.array([Trial]), iTrialLength, PARAMS['Fmin'], PARAMS['Fmax'], PARAMS['frequency'])
    
    if not os.path.isfile(PARAMS['file_path']):
        with open(PARAMS['file_path'], 'w') as out:
            out.write('data_name,oscore\n')
        
    with open(PARAMS['file_path'], 'a+') as out:
            out.write('{},{}\n'.format(PARAMS['data_name'], oscore))
            
    print 'done'
            
            
if __name__ == '__main__':
    main()