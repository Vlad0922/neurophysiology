import numpy as np

# neuron == 0
# interval == 2
def get_spiketrains(nex_data):    
    for v in filter(lambda v: v['Header']['Type'] == 0, nex_data['Variables']):
        yield v['Header']['Name'].lower(), np.array(v['Timestamps'])
        
def get_intervals(nex_data):
    for v in filter(lambda v: v['Header']['Type'] == 2, nex_data['Variables']):
        yield v['Header']['Name'].lower(), list(zip(*v['Intervals']))

        
class Checker(object):
    def __init__(self, merge_intervals=True):
        self.merge_intervals = merge_intervals
    
    def filter_intervals(self, spiketrains, intervals):
        pass
    
    @staticmethod
    def _merge_spiketrains(spiketrain, interval_values):
        st_list = [spiketrain[(spiketrain >= start) & (spiketrain <= end)] for start,end in interval_values]
        offset = st_list[0][~0]

        for i in range(1, len(st_list)):
            st_list[i] = st_list[i][1:] - st_list[i][0] + offset
            offset = st_list[i][~0]

        return np.concatenate(st_list)
    
    def apply_interval(self, spiketrain_name, spiketrain, interval_name, interval_values):
        if self.merge_intervals:
            spikes = self._merge_spiketrains(spiketrain, interval_values)
            yield spiketrain_name, interval_name, spikes
        else:
            for interval in interval_values:
                start, end = interval
                spikes = np.asarray(spiketrain)
                yield spiketrain_name, interval_name, spikes[(spikes >= start) & (spikes <= end)] 
    

class AllfileChecker(Checker):
    def filter_intervals(self, spiketrains, intervals):
        for st_name, st_events in spiketrains:
            yield st_name, 'allfile', st_events

            
class FonChecker(Checker):
    def filter_intervals(self, spiketrains, intervals):
        fon_intervals = [interval for interval in intervals if interval[0] == 'fon']
        
        for st_name, st_events in spiketrains:    
            for interval_name, interval_values in fon_intervals:
                yield from self.apply_interval(st_name, st_events, 'fon', interval_values)
                
class NameChecker(Checker):
    def __init__(self, interval_name, merge_intervals=True):
        super(NameChecker, self).__init__(merge_intervals)
        
        self.interval_name = interval_name
        
    def filter_intervals(self, spiketrains, intervals):
        fon_intervals = [interval for interval in intervals if interval[0] == 'fon']
        
        for st_name, st_events in spiketrains:    
            for interval_name, interval_values, in filter(lambda x: x[0] == self.interval_name, intervals):
                yield from self.apply_interval(st_name, st_events, interval_name, interval_values)


class NeuronChecker(Checker):
    def __init__(self, merge_intervals=True, check_suffix=False):
        super(NeuronChecker, self).__init__(merge_intervals)

        if check_suffix:
            self.name_checker = lambda x,y: y.endswith(x)
        else:
            self.name_checker = lambda x,y: y.startswith(x)
    
    def filter_intervals(self, spiketrains, intervals):
        for st_name, st_events in spiketrains:
            for interval_name, interval_values, in filter(lambda x: self.name_checker(st_name, x[0]), intervals):
                yield from self.apply_interval(st_name, st_events, interval_name, interval_values) 


def apply_intervals(spiketrains, intervals, names_preprocessed=False, fixed_interval_name=None, interval_mode=None):
    interval_names = [interval[0] for interval in intervals]
    
    if not(fixed_interval_name is None):
        checker = NameChecker(fixed_interval_name) 
    elif interval_mode.lower() == 'allfile':
        checker = AllfileChecker()
    else
        checker = FonChecker()
    
    yield from checker.filter_intervals(spiketrains, intervals)