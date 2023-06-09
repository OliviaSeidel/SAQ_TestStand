"""Reconstruct current on a plot using vdd measured from output of integrator slope on the oscilloscope,
corrected capacitance, and a weighted mean of the rtds measured from inputting the known current source onto
each channel"""

import matplotlib.pyplot as plt
import ast
import numpy as np
import statistics

def calc_quartiles_bounds(data):
    q1 = np.percentile(data, 25)
    q2 = np.percentile(data, 50)
    q3 = np.percentile(data, 75)

    iqr = q3 - q1
    upper_bound = q3 + 1.5 * iqr
    lower_bound = q1 - 1.5 * iqr

    return q1, q2, q3, upper_bound, lower_bound
def weighted_average_rtd_bins(rtd_list):
    if len(rtd_list) == 0:
        weightedmean = 0
    else:
        # Define the number of bins
        num_bins = 50

        # Calculate the bin edges
        bin_edges = np.linspace(min(rtd_list), max(rtd_list), num_bins + 1)

        # Bin the values and get the bin counts
        binned_values, _ = np.histogram(rtd_list, bins=bin_edges)

        # Get the bin indices for each value in rtd_list
        bin_indices = np.digitize(rtd_list, bin_edges) - 1

        # Multiply each value by its weight
        weighted_vals = [rtd_list[i] * binned_values[bin_indices[i]] if bin_indices[i] < num_bins else 0 for i in range(len(rtd_list))]

        weighted_sum = sum(weighted_vals)
        sum_weights = sum([binned_values[bin_indices[i]] if bin_indices[i] < num_bins else 0 for i in range(len(rtd_list))])
        weightedmean = weighted_sum/sum_weights

    return weightedmean

# Define the path to the text file
filelist = ["systematics/ch1.txt","systematics/ch2.txt","systematics/ch3.txt","systematics/ch4.txt","systematics/ch5.txt","systematics/ch6.txt","systematics/ch7.txt","systematics/ch8.txt","systematics/ch9.txt","systematics/ch10.txt","systematics/ch11.txt","systematics/ch12.txt","systematics/ch13.txt","systematics/ch14.txt","systematics/ch15.txt","systematics/ch16.txt"]

# Initialize a list to hold the 16 lists of reset time differences
reset_time_diffs = [[] for _ in range(16)]

weightedMeanRtds =[]
avgrtd = []

# Read in the data from the text file and populate the list
toperr=[]
boterr=[]
for filepath in filelist:
    chNum=int(str(filepath).split('ch')[-1].split('.')[0])

    with open(filepath, "r") as f:
        data = ast.literal_eval(f.read())
        for i, rtd_list in enumerate(data):
            if len(rtd_list) > 0:
                reset_time_diffs[i] = rtd_list

    # Loop over each channel and create a histogram plot
    for i, rtd_list in enumerate(reset_time_diffs):

        if len(rtd_list) == 0 or (i+1) != chNum:
            continue

        else:
            q1, q2, q3, upper_bound, lower_bound = calc_quartiles_bounds(rtd_list)
            rtd_list = [x for x in rtd_list if (x <= (upper_bound) and x >= lower_bound)]

            avg_rtd_Weighted = weighted_average_rtd_bins(rtd_list)
            weightedMeanRtds.append(avg_rtd_Weighted)
            std_dev = statistics.stdev(rtd_list)
            one_std_dev_above = avg_rtd_Weighted + std_dev
            one_std_dev_below = avg_rtd_Weighted - std_dev
            toperr.append(one_std_dev_above)
            boterr.append(one_std_dev_below)


VddmVTimes4_1 = [896,896,884,884,884,900,900,888,880,928,972,884,1000,1072,964,932]
VddMeasuredOscilliscope=[241,246,249,247,238,241,247,241,240,248,249,237,256,255,256,250]
CorrectedC = [9.598e-12, 9.643e-12, 1.0045e-11, 1.0045e-11, 1.0045e-11, 1.0089e-11, 9.956e-12, 9.91e-12, 9.636e-12, 1e-11, 9.547e-12, 9.819e-12, 1.004e-11, 1.0112e-11, 1.0083e-11, 9.571e-12]
avgCurrents=[]

toperrCurr=[]
boterrCurr=[]
for i in range (0,len(weightedMeanRtds)):
    rtd = weightedMeanRtds[i]
    c =CorrectedC[i]
    vdd=VddmVTimes4_1[i]*10**(-3)
    avgcurr = ((vdd * c) /rtd)
    avgCurrents.append(avgcurr)
    toperrbar=((vdd * c) /toperr[i])
    boterrbar = ((vdd * c) / boterr[i])
    toperrCurr.append(toperrbar)
    boterrCurr.append(boterrbar)

fig, ax = plt.subplots()
ch =[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
ax.plot(ch,avgCurrents)
ax.axhline(y=1e-9, color='red', linestyle='--', label = '1 nA')
ax.errorbar(ch, avgCurrents, yerr=[boterrCurr, toperrCurr], capsize=4, markersize='2')
ax.set_title("Reconstructed Current")
ax.set_xlabel('Channel Number')
ax.set_ylabel('Current (A)')
ax.legend()
fig.savefig('1nAreconstructedCurrent')

