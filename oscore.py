from __future__ import division

import numpy as np
from numpy.lib.stride_tricks import as_strided

DEFAULT_FREQ = 1./1000

# TODO: REFACTOR IT
# I would have fired the man who wrote it..this is terrible code

def OScoreAC(AutoCorrelogram, LowBoundFrequency, HighBoundFrequency, SamplingFrequency):
    CorrelogramSize = len(AutoCorrelogram)
    HalfCorrelogramSize = int(np.floor(CorrelogramSize/2)+1)

    MirroredAutoCorrelogram = np.flipud(AutoCorrelogram)

    KernelSTD = min(2, 134.0/(1.5*HighBoundFrequency)) * (SamplingFrequency/1000.0)
    GaussianKernel = np.array([1/(np.sqrt(2*np.pi)*KernelSTD)*np.exp(-1.*i*i/(2.0*KernelSTD*KernelSTD))\
                               for i in range(-int(np.round(3*KernelSTD)), int(np.round(3*KernelSTD))+1)])
 
    PaddedAutoCorrelogram =  np.concatenate([MirroredAutoCorrelogram[-len(GaussianKernel):], AutoCorrelogram, MirroredAutoCorrelogram[:len(GaussianKernel)]])
    HighSmoothedAutoCorrelogram = np.convolve(PaddedAutoCorrelogram, GaussianKernel, 'same')
    HighSmoothedAutoCorrelogram = HighSmoothedAutoCorrelogram[len(GaussianKernel):-len(GaussianKernel)]

    KernelSTD =  2 * 134.0/(1.5*LowBoundFrequency) * (SamplingFrequency/1000.0)
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

    LowFrequencyIndex = int(np.round(LowBoundFrequency * SpectrumSize/(SamplingFrequency/2)))
    HighFrequencyIndex = int(np.round(HighBoundFrequency * SpectrumSize/(SamplingFrequency/2)))
    PeakFreqVal = np.max(Spectrum[LowFrequencyIndex:HighFrequencyIndex+1])
    PeakFreqIndex = np.where(Spectrum == PeakFreqVal)[0][0]

    if(SpectrumIntegral > 0):
        OscScore = PeakFreqVal / SpectrumIntegral;
        OscFreq = ((PeakFreqIndex-1+LowFrequencyIndex-1) * (SamplingFrequency/2))/SpectrumSize
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
    return T.dot(p[maxlag:].conj())


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
    return T.dot(px)


def oscore_spikes(SpikeTimesPerTrials, TrialLength, LowBoundFrequency, HighBoundFrequency, SamplingFrequency):
    CorrelationWindow = int(np.floor(np.power(2, np.ceil(max(np.log2(3*SamplingFrequency/LowBoundFrequency),np.log2(SamplingFrequency/4.0))))))
    AutoCorrelogramSize = 2 * CorrelationWindow + 1

    TrialCount = len(SpikeTimesPerTrials)
    TrialAutoCorrelograms = np.zeros((TrialCount,AutoCorrelogramSize))
    TrialOScores = np.zeros(TrialCount)

    for i in range(TrialCount):
        Trial = np.zeros(TrialLength)
        Trial[SpikeTimesPerTrials[i]-1] = 1
        TrialAutoCorrelograms[i,:] = crosscorrelation(Trial, Trial, CorrelationWindow)
        TrialOScores[i] = OScoreAC(TrialAutoCorrelograms[i,:], LowBoundFrequency, HighBoundFrequency, SamplingFrequency)

    AutoCorrelogram = np.sum(TrialAutoCorrelograms, axis=0)

    OSc = OScoreAC(AutoCorrelogram, LowBoundFrequency, HighBoundFrequency, SamplingFrequency)
    return OSc
