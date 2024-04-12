import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

# calculate the times and magnitudes of spikes from an SSR trace
def toSpikes(trace, sampleRate):
    search_range = int(0.0016 * sampleRate)
    maxes = sp.signal.find_peaks(trace)[0]  # find the time locations of the relative maximums in the trace
    # inner function to calculate one spike from a maximum found
    def calcSpike(x1):
        x2 = x1 + np.argmin(trace[x1: x1 + search_range])  # calculate position of minimum value
        x1 = x2 - search_range + np.argmax(trace[x2 - search_range : x2])   # recalculate position of maximum by checking if there were any bigger maximums nearby
        avg = (x1 + x2) / 2     # set position of spike in the middle of the min and max
        spike = int(trace[x1] - trace[x2])  # calculate magnitude of spike as difference between max and min amplitudes
        if spike > 0:   # if max is correctly higher than the min, return a row as the spike's time and magnitude
            return [avg / sampleRate, spike]
        else:   # else do not return a spike
            return None
    # !!! not unique, can be duplicates? if multiple local maximums near each other, consider using np.unqiue(, axis=0)
        # fits AutoSpike histogram better?
    # calculate all spikes based on each relative maximum and return in a 2D array with the first row as the time locations and the second row as the magnitudes
    ###return np.array([s for s in [calcSpike(x1) for x1 in maxes if x1 + SEARCH_RANGE < len(trace) and x1 > SEARCH_RANGE] if s != None]).T
    spikes = np.unique(np.array([s for s in [calcSpike(x1) for x1 in maxes if x1 + search_range < len(trace) and x1 > search_range] if s != None]), axis=0).T
    spike_dict = {}
    last = -1
    for i in range(len(spikes[0])):     # removing duplicate spikes if there are spikes less than a milisecond apart
        time = spikes[0][i]
        if abs(last - time) < 0.001:
            spike_dict[last].append(spikes[1][i])
        else:
            last = time
            spike_dict[last] = [spikes[1][i]]
    result = []
    for t in spike_dict.keys():
        if len(spike_dict[t]) > 1:
            result.append([t, max(spike_dict[t])])
        else:
            result.append([t, spike_dict[t][0]])
    return np.array(result).T

def spikesInRange(spikes, start, end, min_height, max_height):
    def check_spike(i):    
        time = spikes[0][i]
        if time >= start and time <= end:
            height = spikes[1][i]
            if height >= min_height and height <= max_height:
                return [time, height]
        return None
    return np.array([s for s in [check_spike(i) for i in range(len(spikes[0]))] if s != None]).T    # use a faster search method

def plotSpikes(spikes, color):
    if spikes.size > 0:
        x = spikes[0]
        y = spikes[1]
        y_min = -y/2
        y_max = -y_min
        plt.vlines(x, y_min, y_max, color=color)

def plotVlines(xs, color, linestyle='solid'):
    for x in xs:
        plt.axvline(x, color=color, linestyle=linestyle)

