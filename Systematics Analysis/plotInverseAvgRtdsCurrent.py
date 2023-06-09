"""Based on 3 different input current into each channel, this outputs a plot with 16 subplots with the average rtd
per input current. Done using box/whisker outlier removal method, and weighted mean"""

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
filelist1 = ["txtfiles1nA/ch1.txt","txtfiles1nA/ch2.txt","txtfiles1nA/ch3.txt","txtfiles1nA/ch4.txt","txtfiles1nA/ch5.txt","txtfiles1nA/ch6.txt","txtfiles1nA/ch7.txt","txtfiles1nA/ch8.txt","txtfiles1nA/ch9.txt","txtfiles1nA/ch10.txt","txtfiles1nA/ch11.txt","txtfiles1nA/ch12.txt","txtfiles1nA/ch13.txt","txtfiles1nA/ch14.txt","txtfiles1nA/ch15.txt","txtfiles1nA/ch16.txt"]
filelist2 = ["txtfiles2nA/ch1.txt","txtfiles2nA/ch2.txt","txtfiles2nA/ch3.txt","txtfiles2nA/ch4.txt","txtfiles2nA/ch5.txt","txtfiles2nA/ch6.txt","txtfiles2nA/ch7.txt","txtfiles2nA/ch8.txt","txtfiles2nA/ch9.txt","txtfiles2nA/ch10.txt","txtfiles2nA/ch11.txt","txtfiles2nA/ch12.txt","txtfiles2nA/ch13.txt","txtfiles2nA/ch14.txt","txtfiles2nA/ch15.txt","txtfiles2nA/ch16.txt"]
filelist3 = ["txtfiles3nA/ch1.txt","txtfiles3nA/ch2.txt","txtfiles3nA/ch3.txt","txtfiles3nA/ch4.txt","txtfiles3nA/ch5.txt","txtfiles3nA/ch6.txt","txtfiles3nA/ch7.txt","txtfiles3nA/ch8.txt","txtfiles3nA/ch9.txt","txtfiles3nA/ch10.txt","txtfiles3nA/ch11.txt","txtfiles3nA/ch12.txt","txtfiles3nA/ch13.txt","txtfiles3nA/ch14.txt","txtfiles3nA/ch15.txt","txtfiles3nA/ch16.txt"]
currentFiles= ['drive-download-20230604T191134Z-001/1nA/currentSourceUncertainty/1nA/current_time0.txt','drive-download-20230604T191134Z-001/2nA/currentSourceUncertainty/current_time0.txt','drive-download-20230604T191134Z-001/3nA/currentSourceUncertainty/current_time0.txt']

currents=[]
std_devs=[]
for filepath in currentFiles:
    values_list = []
    with open(filepath, "r") as f:
        for line in f:
            # Split the line by space and take the first column
            columns = line.split()
            if columns:
                value = columns[0]
                values_list.append(float(value))
    std_devs.append(3*np.std(values_list))
    currents.append(values_list)
print(std_devs)

# Initialize a list to hold the 16 lists of reset time differences
reset_time_diffs = [[] for _ in range(16)]

avgRtd1nA = []
std_dev_1nA=0
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
            rtd_list = [1/rtd for rtd in rtd_list]
            avgRtd1nA.append(weighted_average_rtd_bins(rtd_list))
            std_dev_1nA= 2*statistics.stdev(rtd_list)


avgRtd2nA = []
std_dev_2nA = 0
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
            rtd_list = [1 / rtd for rtd in rtd_list]

            avgRtd2nA.append(weighted_average_rtd_bins(rtd_list))
            std_dev_2nA = 2*statistics.stdev(rtd_list)


avgRtd3nA = []
std_dev_3nA =0
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
            rtd_list = [1 / rtd for rtd in rtd_list]
            avgRtd3nA.append(weighted_average_rtd_bins(rtd_list))
            std_dev_3nA = 2*statistics.stdev(rtd_list)

fig, axs = plt.subplots(4, 4, figsize=(12, 12))
axs = axs.ravel()
for i in range(0,16):
    xnumRtds=[round(avgRtd1nA[i],4),round(avgRtd2nA[i],4),round(avgRtd3nA[i],4)]
    xer= [[std_dev_1nA,std_dev_2nA,std_dev_3nA],[std_dev_1nA,std_dev_2nA,std_dev_3nA] ]
    y=[1e-9,2e-9,3e-9]
    # need to plot the error not the actual value (point-error)...matplotlib is smart
    yer = [[std_devs[0],std_devs[1], std_devs[2]], [std_devs[0],std_devs[1], std_devs[2]]]
    axs[i].plot(xnumRtds, y)
    axs[i].errorbar(xnumRtds,y, yerr=yer,xerr=xer, fmt='o', markersize=2, linewidth=0.01, capsize=8, color='black', ecolor='red')
    axs[i].set_xlabel('1/(Avg Rtd) (1/s)')
    axs[i].set_ylabel('Input I (A)')
    axs[i].set_yticks(y)
    axs[i].set_xticks(xnumRtds)
    axs[i].set_title('Ch ' + str(i) + '1/(RTD) per Input I')
    axs[i].figure.show()
fig.tight_layout()
fig.savefig('currentVsRtdsAllCh'+'.png', dpi=500)

