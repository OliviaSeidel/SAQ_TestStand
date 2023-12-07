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

import  os
filepaths_list = []

# Get the current directory where this script is located
current_dir = os.path.dirname(os.path.abspath(__file__))

filepaths_list = []

filepaths_list = [os.path.join(current_dir, file) for file in os.listdir(current_dir) if file.endswith(".txt")]

# Sort the list of file paths alphabetically
filepaths_list = sorted(filepaths_list)

# Define the path to the text file
filelist1 = []
filelist2 = []
filelist3 =[]
filelist4 =[]

for file in filepaths_list:
    if "0p01nA" in file:
        filelist1.append(file)
    if "0p07nA" in file:
        filelist2.append(file)
    if "0p1nA" in file:
        filelist3.append(file)
    if "0p5nA" in file:
        filelist4.append(file)
#0.01,0.07,0.1,0.5nA
filelists=[filelist1,filelist2,filelist3,filelist4]
datalists=[]
averageRtdLists=[]
n=0
for filelist in filelists: #go through each current
    thisCurrentsData = [[] for _ in range(14)]  #for each current there are 14 channels
    thisCurrentAverages = [[] for _ in range(14)]
    for thisfilepath in filelist: #go t0 each channel
        with open(thisfilepath, "r") as f:
            data = ast.literal_eval(f.read()) #get the data which has 16 channels but you only need 1
            chnum = int(thisfilepath.split("Ch")[1].split("_")[0])-1 #get the channel number

            times=[]
            for timeStamp in data[chnum]: #convert each one to new time
                times.append(timeStamp * 200 / 30e6)
            rtd_list = []
            for index in range(9, len(times)): #skip the 0 packets in front
                rtd = times[index] - times[index - 1]
                rtd_list.append(rtd)
            if not rtd_list:
                print("The list is empty.")
            thisCurrentAverages[chnum]=np.average(rtd_list)
            thisCurrentsData[chnum]=rtd_list
    datalists.append(thisCurrentsData)
    averageRtdLists.append(thisCurrentAverages)
    n+=1


fig, axs = plt.subplots(4, 4, figsize=(12, 12))
axs = axs.ravel()

for i in range(14):
    x = [0.01, 0.07, 0.1, 0.5]
    y = [1 / averageRtdLists[0][i], 1 / averageRtdLists[1][i], 1 / averageRtdLists[2][i], 1 / averageRtdLists[3][i]]

    # Perform linear regression
    coefficients = np.polyfit(x, y, 1) #returns slope and y intercept of fit as coefficients[0] and coefficients[1]
    print(1/coefficients[0])
    polynomial = np.poly1d(coefficients)
    line_of_best_fit = polynomial(x)

    # Plot the data points
    axs[i].plot(x, y, 'o', label=f'ChargePerReset: {1/coefficients[0]:.4f}')
    axs[i].plot(x, y, 'o', label='Data points')

    # Plot the line of best fit
    axs[i].plot(x, line_of_best_fit, label=f'Line of best fit: {coefficients[0]:.4f}x + {coefficients[1]:.4f} ')


    # Set labels and title for the subplot
    axs[i].set_ylabel('1/(Avg Rtd) (1/s)')
    axs[i].set_xlabel('Input I (A)')
    axs[i].set_title('Ch ' + str(i) + '  1/(RTD) per Input I')

    # Show legend with smaller font size
    axs[i].legend(fontsize=6)

# Adjust layout for better spacing
plt.tight_layout()

fig.savefig('Calibration11152023'+'.png', dpi=500)

