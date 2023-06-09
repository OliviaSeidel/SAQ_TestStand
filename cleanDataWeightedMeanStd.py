
"""This takes an input of 16 different .txt files, each containing resets on all 16 channels. It plots each consecutive channel.
This is used to plot the calibration results. Inputting current on one channel, and recording the resets on all channels
consecutively.

It outputs a plot with 16 subplots with cleaned data (using box/whisker outlier removal method) and plots the wighted mean
and 1 standard deviation from the mean"""

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
        weighted_mean = 0
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
        weighted_mean = weighted_sum / sum_weights

    return weighted_mean

# Define the path to the text file
filelist = ["txtfiles3nA/ch{}.txt".format(i) for i in range(1, 17)]

reset_time_diffs = [[] for _ in range(16)]

avgrtd = []
fig, axs = plt.subplots(4, 4, figsize=(12, 12))
axs = axs.ravel()

for filepath in filelist:
    chNum = int(filepath.split('ch')[-1].split('.')[0])

    with open(filepath, "r") as f:
        data = ast.literal_eval(f.read())
        for i, rtd_list in enumerate(data):
            if len(rtd_list) > 0:
                reset_time_diffs[i] = rtd_list

    for i, rtd_list in enumerate(reset_time_diffs):
        if len(rtd_list) == 0 or (i + 1) != chNum:
            continue
        else:
            q1, q2, q3, upper_bound, lower_bound = calc_quartiles_bounds(rtd_list)
            rtd_list = [x for x in rtd_list if (x <= upper_bound and x >= lower_bound)]

            avg_rtd_weighted = weighted_average_rtd_bins(rtd_list)
            std_dev = statistics.stdev(rtd_list)
            one_std_dev_above = avg_rtd_weighted + std_dev
            one_std_dev_below = avg_rtd_weighted - std_dev

            axs[i].hist(rtd_list, bins=50)
            axs[i].axvline(avg_rtd_weighted, color='red', linestyle='--', label='weighted mean')
            axs[i].axvline(one_std_dev_above, color='blue', linestyle='--', label='1 std above')
            axs[i].axvline(one_std_dev_below, color='green', linestyle='--', label='1 std below')
            axs[i].set_xlabel("rtd")
            axs[i].set_ylabel("Frequency")
            axs[i].legend(fontsize='6')
            axs[i].set_title("Ch{} (Avg rtd: {:.2e})".format(i + 1, avg_rtd_weighted))

fig.tight_layout()
fig.savefig('16rtd')
