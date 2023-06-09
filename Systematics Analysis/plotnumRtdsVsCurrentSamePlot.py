"""Based on 3 diffierent input current into each channel, this outputs a plot with 16 subplots with the number of rtds
per input current"""

import matplotlib.pyplot as plt
import ast
import numpy as np

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
filelist1 = ["txtfiles1nA/ch1.txt","txtfiles1nA/ch2.txt","txtfiles1nA/ch3.txt","txtfiles1nA/ch4.txt","txtfiles1nA/ch5.txt","txtfiles1nA/ch6.txt","txtfiles1nA/ch7.txt","txtfiles1nA/ch8.txt","txtfiles1nA/ch9.txt","txtfiles1nA/ch10.txt","txtfiles1nA/ch11.txt","txtfiles1nA/ch12.txt","txtfiles1nA/ch13.txt","txtfiles1nA/ch14.txt","txtfiles1nA/ch15.txt","txtfiles1nA/ch16.txt"]
filelist2 = ["txtfiles2nA/ch1.txt","txtfiles2nA/ch2.txt","txtfiles2nA/ch3.txt","txtfiles2nA/ch4.txt","txtfiles2nA/ch5.txt","txtfiles2nA/ch6.txt","txtfiles2nA/ch7.txt","txtfiles2nA/ch8.txt","txtfiles2nA/ch9.txt","txtfiles2nA/ch10.txt","txtfiles2nA/ch11.txt","txtfiles2nA/ch12.txt","txtfiles2nA/ch13.txt","txtfiles2nA/ch14.txt","txtfiles2nA/ch15.txt","txtfiles2nA/ch16.txt"]
filelist3 = ["txtfiles3nA/ch1.txt","txtfiles3nA/ch2.txt","txtfiles3nA/ch3.txt","txtfiles3nA/ch4.txt","txtfiles3nA/ch5.txt","txtfiles3nA/ch6.txt","txtfiles3nA/ch7.txt","txtfiles3nA/ch8.txt","txtfiles3nA/ch9.txt","txtfiles3nA/ch10.txt","txtfiles3nA/ch11.txt","txtfiles3nA/ch12.txt","txtfiles3nA/ch13.txt","txtfiles3nA/ch14.txt","txtfiles3nA/ch15.txt","txtfiles3nA/ch16.txt"]


# t = lambda x: x * SAQ_DIV / ZYBO_FRQ
# Initialize a list to hold the 16 lists of reset time differences
reset_time_diffs = [[] for _ in range(16)]

numRtds1nA = []
for filepath in filelist1:
    chNum=int(str(filepath).split('ch')[-1].split('.')[0])
    with open(filepath, "r") as f:
        data = ast.literal_eval(f.read())
        for i, rtd_list in enumerate(data):
            if len(rtd_list) > 0:
                reset_time_diffs[i] = rtd_list
    for i, rtd_list in enumerate(reset_time_diffs):
        if len(rtd_list) == 0 or (i+1) != chNum:
            continue
        else:
            q1, q2, q3, upper_bound, lower_bound = calc_quartiles_bounds(rtd_list)
            rtd_list = [x for x in rtd_list if (x <= (upper_bound) and x >= lower_bound)]
            numRtds1nA.append(len(rtd_list))

numRtds2nA = []
for filepath in filelist2:
    chNum=int(str(filepath).split('ch')[-1].split('.')[0])
    with open(filepath, "r") as f:
        data = ast.literal_eval(f.read())
        for i, rtd_list in enumerate(data):
            if len(rtd_list) > 0:
                reset_time_diffs[i] = rtd_list
    for i, rtd_list in enumerate(reset_time_diffs):
        if len(rtd_list) == 0 or (i+1) != chNum:
            continue
        else:
            q1, q2, q3, upper_bound, lower_bound = calc_quartiles_bounds(rtd_list)
            rtd_list = [x for x in rtd_list if (x <= (upper_bound) and x >= lower_bound)]
            numRtds2nA.append(len(rtd_list))

numRtds3nA = []
for filepath in filelist3:
    chNum=int(str(filepath).split('ch')[-1].split('.')[0])
    with open(filepath, "r") as f:
        data = ast.literal_eval(f.read())
        for i, rtd_list in enumerate(data):
            if len(rtd_list) > 0:
                reset_time_diffs[i] = rtd_list
    for i, rtd_list in enumerate(reset_time_diffs):
        if len(rtd_list) == 0 or (i+1) != chNum:
            continue
        else:
            q1, q2, q3, upper_bound, lower_bound = calc_quartiles_bounds(rtd_list)
            rtd_list = [x for x in rtd_list if (x <= (upper_bound) and x >= lower_bound)]
            numRtds3nA.append(len(rtd_list))
fig, axs = plt.subplots(4, 4, figsize=(12, 12))
axs = axs.ravel()
for i in range(0,16):
    xnumRtds=[numRtds1nA[i],numRtds2nA[i],numRtds3nA[i]]
    y=[1,2,3]
    axs[i].plot(xnumRtds,y,marker = 'o')
    axs[i].set_xlabel('Number of Rtds')
    axs[i].set_ylabel('Input I (nA)')
    axs[i].set_yticks(y)
    axs[i].set_xticks(xnumRtds)
    axs[i].set_title('Ch ' + str(i) + ' Rtds per Input I')
fig.tight_layout()
fig.savefig('currentVsRtdsAllCh'+'.png')

