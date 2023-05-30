import numpy as np
import copy

from scipy.io import wavfile
from scipy.signal import detrend, resample, butter, filtfilt

def ms_to_samples(ms, fs):
    """
    Function for converting a given number of milliseconds (and a sample rate `fs`) to samples.

    >>> ms_to_samples(107, 44100)
    4719
    """
    return int(np.floor(ms / 1000 * fs) + 1)

def samples_to_ms(sample, fs):
    """
    Function for converting a given number of samples (and a sample rate `fs`) to milliseconds. 
    """
    return sample * 1000 / fs

def apply_fades(sig, ms, fs):
    """
    Function for applying an exponential fade to a given signal.
    """
    fade_in = np.power(np.linspace(0, 1, 1 + int(np.round(ms / 1000 * fs))), 2)
    fade_out = fade_in[::-1]
    sig[:len(fade_in)] = sig[:len(fade_in)] * fade_in
    sig[-len(fade_out):] = sig[-len(fade_out):] * fade_out

    return sig

def se_resample(samples, source_rate, target_rate):
	resample_factor = float(target_rate) / float(source_rate)
	resampled=resample(samples, int(len(samples) * resample_factor))
	tt = np.linspace(0, (len(resampled)) / (1.0 * target_rate), len(resampled))
	
	return resampled, tt

def freq2nyq(freq, fs):
	# frequencies for filter design must be expressed as a fraction of Nyquist
	# nyquist is half the sample rate
	nyquist = fs / 2.0
	return freq / nyquist

def nyq2freq(nyq, fs):
	nyquist = fs / 2.0
	return nyq * nyquist

def build_filter_b(Wn):
	b, a = butter(2, Wn, btype='band', analog=False, output='ba')
	return b, a

def band_pass_sharp(signal, fs, lower, upper):
	Wn1 = [freq2nyq(lower, fs), freq2nyq(upper, fs)]
	b1, a1 = build_filter_b(Wn1)
	filtered_signal1 = filtfilt(b1, a1, signal)
	filtered_signal2 = filtfilt(b1, a1, filtered_signal1)
	return filtered_signal2

def save_samples_to_file(samples, filename, fs):
    wavfile.write(filename, rate=fs, data=samples.astype(np.float32))

def midi2freq(midi_number):
    return (440 / 32) * (2**((midi_number - 9) / 12))

def freq2midi(f0):
	if f0 <= 0:
		return 0
	else:
		exp = result=12 * (np.log2(f0) - np.log2(440)) + 69
		if exp < 0:
			return 0 # get's weird log2 values otherwise...
		else:
			return exp

def diff(x): 
    return [j - i for i, j in zip(x[:-1], x[1:])]

def diff_across(x, y):
	return [j - i for i, j in zip(x, y)]

def to_contour(vals):
	seg_vals = copy.copy(vals)
	seg_vals = [round(x) for x in seg_vals]
	value_dict = dict()

	for i, this_val in zip(range(0, len(sorted(set(seg_vals)))), sorted(set(seg_vals))):
		value_dict[this_val] = str(i)

	for i, this_val in enumerate(seg_vals):
		for this_key in value_dict:
			if this_val == this_key:
				seg_vals[i] = value_dict[this_key]

	return [int(val) for val in seg_vals]

def uds(x):
	if x>0.5:
		return 1
	elif x<-0.5:
		return -1
	else:
		return 0

def uds_contour_interval(intervals):
	return [uds(x) for x in intervals]

def save_samples_to_file(samples, filename, fs):
    wavfile.write(filename, rate=fs, data=np.array(samples, dtype=np.float32))

# helper function to evaluate if two segments (a,b) and (x,y) overlap
# 0 it segements (a,b) and (x,y) don't overlap
def overlap_segments_score(a,b,x,y):
    if ((b<a) or (y<x)):
        return -1 # this is not a sgement
    u=max(a,x)
    v=min(b,y)
    score = (v-u)/min(b-a,y-x)
    if ((x==b) or (a==y)):
        return 0.5 # just touching!
    if score<=0:
        return 0 # not overlapping
    else:
        return 1 # overlapping
