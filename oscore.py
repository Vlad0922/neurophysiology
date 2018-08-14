from __future__ import division

import numpy as np
from numpy.lib.stride_tricks import as_strided

DEFAULT_FREQ = 1./1000

# TODO: REFACTOR IT
# I would have fired the man who wrote it..this is terrible code

def OScoreAC(AutoCorrelogram, low_bound, high_bound, freq):
    CorrelogramSize = len(AutoCorrelogram)
    HalfCorrelogramSize = int(np.floor(CorrelogramSize/2)+1)

    MirroredAutoCorrelogram = np.flipud(AutoCorrelogram)

    KernelSTD = min(2, 134.0/(1.5*high_bound)) * (freq/1000.0)
    GaussianKernel = np.array([1/(np.sqrt(2*np.pi)*KernelSTD)*np.exp(-1.*i*i/(2.0*KernelSTD*KernelSTD))\
                               for i in range(-int(np.round(3*KernelSTD)), int(np.round(3*KernelSTD))+1)])
 
    PaddedAutoCorrelogram =  np.concatenate([MirroredAutoCorrelogram[-len(GaussianKernel):], AutoCorrelogram, MirroredAutoCorrelogram[:len(GaussianKernel)]])
    HighSmoothedAutoCorrelogram = np.convolve(PaddedAutoCorrelogram, GaussianKernel, 'same')
    HighSmoothedAutoCorrelogram = HighSmoothedAutoCorrelogram[len(GaussianKernel):-len(GaussianKernel)]

    KernelSTD =  2 * 134.0/(1.5*low_bound) * (freq/1000.0)
    GaussianKernel = np.array([1/(np.sqrt(2*np.pi)*KernelSTD)*np.exp(-1.*i*i/(2.0*KernelSTD*KernelSTD))\
                               for i in range(-int(np.round(3*KernelSTD)), int(np.round(3*KernelSTD))+1)])

    PaddedAutoCorrelogram =  np.concatenate([MirroredAutoCorrelogram[-len(GaussianKernel):], AutoCorrelogram, MirroredAutoCorrelogram[:len(GaussianKernel)]])
    LowSmoothedAutoCorrelogram = np.convolve(PaddedAutoCorrelogram, GaussianKernel, 'same')
    LowSmoothedAutoCorrelogram = LowSmoothedAutoCorrelogram[len(GaussianKernel):-len(GaussianKernel)]

    Min = min(np.min(LowSmoothedAutoCorrelogram), 0)
    Max = np.max(LowSmoothedAutoCorrelogram)
    if (Max - Min) != 0:
        Ratio = CorrelogramSize/(Max - Min)
    else:
        if (Max == 0):
            return -1.
        else:
            Ratio = 1

    LimitSlopeLeft = np.tan(10.0 * np.pi / 180.0);
    Derivative = np.ediff1d(LowSmoothedAutoCorrelogram)*Ratio

    LeftCutIndex = np.where(Derivative[:HalfCorrelogramSize-1] < LimitSlopeLeft)[0][-1]
    RightCutIndex = HalfCorrelogramSize + (HalfCorrelogramSize-LeftCutIndex - 2)

    AutoCorrelogramWithoutPeak = HighSmoothedAutoCorrelogram
    AutoCorrelogramWithoutPeak[LeftCutIndex:RightCutIndex+1] = np.mean([HighSmoothedAutoCorrelogram[LeftCutIndex], HighSmoothedAutoCorrelogram[RightCutIndex]])

    BlackManWin = np.blackman(len(AutoCorrelogramWithoutPeak))

    Spectrum = np.fft.fft(AutoCorrelogramWithoutPeak*BlackManWin, int(2**np.ceil(np.log2(HalfCorrelogramSize))))

    Spectrum = np.abs(Spectrum)
    SpectrumSize = len(Spectrum)//2
    Spectrum = Spectrum[:SpectrumSize]
    SpectrumIntegral = np.mean(Spectrum)

    LowFrequencyIndex = int(np.round(low_bound * SpectrumSize/(freq/2)))
    HighFrequencyIndex = int(np.round(high_bound * SpectrumSize/(freq/2)))
    PeakFreqVal = np.max(Spectrum[LowFrequencyIndex:HighFrequencyIndex+1])
    PeakFreqIndex = np.where(Spectrum == PeakFreqVal)[0][0]

    if(SpectrumIntegral > 0):
        OscScore = PeakFreqVal / SpectrumIntegral;
        OscFreq = ((PeakFreqIndex-1+LowFrequencyIndex-1) * (freq/2))/SpectrumSize
    else:
        OscScore = -1; 
    
    return OscScore


def _check_arg(x, xname):
    """
    Just a simple dimension check for array
    """
    x = np.asarray(x)
    if x.ndim != 1:
        raise ValueError('%s must be one-dimensional.' % xname)
    return x

def chunking_dot(big_matrix, small_matrix, chunk_size=100):
    # Make a copy if the array is not already contiguous
    small_matrix = np.ascontiguousarray(small_matrix)
    is_flat = (len(small_matrix.shape) == 1)
    if is_flat:        
        small_matrix = small_matrix.reshape((-1,1))

    R = np.empty((big_matrix.shape[0], small_matrix.shape[1]))
    for i in range(0, R.shape[0], chunk_size):
        end = i + chunk_size
        R[i:end] = np.dot(big_matrix[i:end], small_matrix)

    if is_flat:
        R = R.ravel()

    return R

def autocorrelation(x, maxlag):
    """
    Autocorrelation with a maximum number of lags.

    `x` must be a one-dimensional numpy array.

    This computes the same result as
        numpy.correlate(x, x, mode='full')[len(x)-1:len(x)+maxlag]

    The return value has length maxlag + 1.
    """
    x = _check_arg(x, 'x')
    p = np.pad(x.conj(), maxlag, mode='constant')
    T = as_strided(p[maxlag:], shape=(maxlag+1, len(x) + maxlag),
                   strides=(-p.strides[0], p.strides[0]))
    return chunking_dot(T, p[maxlag:].conj())


def crosscorrelation(x, y, maxlag):
    """
    Cross correlation with a maximum number of lags.

    `x` and `y` must be one-dimensional numpy arrays with the same length.

    This computes the same result as
        numpy.correlate(x, y, mode='full')[len(a)-maxlag-1:len(a)+maxlag]

    The return vaue has length 2*maxlag + 1.
    """
    x = _check_arg(x, 'x')
    y = _check_arg(y, 'y')
    py = np.pad(y.conj(), 2*maxlag, mode='constant')
    T = as_strided(py[2*maxlag:], shape=(2*maxlag+1, len(y) + 2*maxlag),
                   strides=(-py.strides[0], py.strides[0]))
    px = np.pad(x, maxlag, mode='constant')
    return chunking_dot(T, px)


def oscore_spikes(spikes_per_trials, trial_len, low_bound, high_bound, freq):
    corr_window = int(np.floor(np.power(2, np.ceil(max(np.log2(3*freq/low_bound), np.log2(freq/4.0))))))
    autocorr_size = 2 * corr_window + 1

    trial_count = len(spikes_per_trials)
    auto_cors = np.zeros((trial_count, autocorr_size))
    triels_oscores = np.zeros(trial_count)

    for i in range(trial_count):
        trial = np.zeros(trial_len)
        trial[spikes_per_trials[i]-1] = 1
        auto_cors[i,:] = crosscorrelation(trial, trial, corr_window)
        triels_oscores[i] = OScoreAC(auto_cors[i,:], low_bound, high_bound, freq)

    autocor = np.sum(auto_cors, axis=0)

    return OScoreAC(autocor, low_bound, high_bound, freq)
