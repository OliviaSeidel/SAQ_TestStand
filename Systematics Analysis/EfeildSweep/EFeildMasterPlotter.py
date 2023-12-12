import matplotlib.pyplot as plt
import ast
import numpy as np
import  os
def sum_list_values(lst):
    total = 0
    for value in lst:
        total += value
    return total

def plotDt(data_dict, errsAbove, errsBelow, fiducialize_num_start=0, fiducialize_num_end=16,DtUnits='cm^s/s', area_normalize=True, mean_from_data=False, fixed_mean = False):
    """Input a dictionary with pressure/field keys
    Output plots of Dt, sigma, fitted gaussians peak normalized, and fitted gaussians non peak normalized"""

    # Initialize variables
    plt.figure()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'lime', 'pink', 'silver', 'gray', 'maroon', 'olive', 'teal', 'navy']
    standardDevs = []
    pressures = []
    plotting_sets = []
    i = 0

    mm = [0.54, 1.3, 2, 2.8, 3.9, 4.9, 5.9, 7.4, 9.9, 12.4, 14.9, 19.9, 24.85, 29.9, 39.8, 49.9]
    mm = mm[:fiducialize_num_end]
    mm = mm[fiducialize_num_start:]

    # Get the charge values for each pressure into a list
    for key, value in data_dict.items():
        plotting_data = None

        pressure = key
        pressure = round(float(pressure) * 51.715, 1)
        pressures.append(pressure)

        if area_normalize:
            plotting_data, area_err = areaNormalize(value, fiducialize_num_start)

        if not area_normalize:
            plotting_data = value

        plotting_sets.append(plotting_data)
        i += 1

    # If you want to use the method where you extract the mean from data:
    if mean_from_data:
        totChargePerPress, standardDevs, mean, numerator = calculate_standard_deviations_mean(mm, plotting_sets, meanVal = 'DataSet')
        plot_sigma(pressures, standardDevs, numerator, totChargePerPress, errsAbove, errsBelow, meanFrom= 'DataSet')

    # If you want to use a fixed mean where you specify how far off from the center the spot center is in mm
    if not mean_from_data:

        # From the definiton call, grab the specified spot center
        mean = fixed_mean

        # Calculate the standard deviations
        totChargePerPress, standardDevs, mean, numerator = calculate_standard_deviations_mean(mm, plotting_sets, meanVal = mean)
        plot_sigma(pressures, standardDevs,  numerator, totChargePerPress, errsAbove, errsBelow, meanFrom = mean)


    return standardDevs

def calculate_standard_deviations_mean(mm, plotting_sets, meanVal= 'DataSet'):

    standardDevs = []

    for i in range(0, len(plotting_sets)):
        # Convert the total-charge per channel data for the current pressure to a NumPy array
        ydata_norm = np.array(plotting_sets[i])

        # Define the number of bins for the histogram to get sigma from
        bins = 300

        # Calculate the histogram with specified bins and use 'ydata_norm' as weights
        counts, bin_edges = np.histogram(mm, bins=bins, weights=ydata_norm)

        # Filter out bins with zero counts
        non_zero_bins = counts > 0

        # Calculate the bin centers for non-zero bins
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_centers = bin_centers[non_zero_bins]

        # Calculate the weighted sum of bin_centers multiplied by counts
        weighted_sum = np.sum(bin_centers * counts[non_zero_bins])

        # Calculate the total sum of counts
        total_count = np.sum(counts)

        # Calculate the mean based on if you have a fixed mean or are getting it from the dataset
        if meanVal == 'DataSet':
            mean = weighted_sum / total_count
            print("Mean of histogrammed data: ", mean)

        else:
            mean = meanVal
            print("Set Mean: ", mean)

        # The following lines find the standard deviation from the histogram
        # Can find this information here: https://www.statology.org/histogram-standard-deviation/

        # Calculate the squared difference between each bin center and the mean
        centerminusmean = ((bin_centers - mean) ** 2)

        # Multiply each squared difference by the corresponding counts for non-zero bins
        countstimesdiff = counts[non_zero_bins] * centerminusmean

        # Calculate the total charge per press by summing the counts for non-zero bins
        totChargePerPress = sum_list_values(counts[non_zero_bins])

        # Define nMinus1 as the total charge per press
        nMinus1 = totChargePerPress

        # Calculate the numerator by summing the products of counts times squared differences
        numerator = sum_list_values(countstimesdiff)

        # Calculate the standard deviation of 'mm' using the formula for a sample
        std_mm = np.sqrt(numerator / nMinus1)

        # Append the standard deviation for this pressure to your total standard deviations
        standardDevs.append(std_mm)

    # Return components for error propagation, like the total Q across all channels seen per pressure, numerator, mean
    return totChargePerPress, standardDevs,mean, numerator

def plot_sigma(pressures, standardDevs,  numerator, totChargePerPress, errsAbove, errsBelow, meanFrom= 'DataSet'):

    YerrsAbove = []
    YerrsBelow = []

    # Go through each of the standard deviations, propagate the error
    for i, std_mm in enumerate(standardDevs):

        # Done using error propagation rules. http://science.clemson.edu/physics/labs/tutorials/errorp/index.html
        # The error in sigma is based on the total charge, so this is propagated through sigma calculation
        # Checked using https://www.wolframalpha.com/widgets/view.jsp?id=8ac60957610e1ee4894b2cd58e753
        errAbove = (numerator * np.sum(errsAbove[i])) / (20 * np.sqrt(numerator/totChargePerPress) * (totChargePerPress ** 2))
        errBelow = (numerator * np.sum(errsBelow[i])) / (20 * np.sqrt(numerator/totChargePerPress) * (totChargePerPress ** 2))

        YerrsAbove.append(errAbove)
        YerrsBelow.append(errBelow)

    # Divide by sqrt of drift length to get the standardized units
    standardDevs = [x/np.sqrt(100) for x in standardDevs] #100mm is 10 cm for 10 cm drift
    plt.figure()
    #plt.errorbar(pressures, standardDevs, yerr=(YerrsBelow, YerrsAbove), fmt='none', capsize=6, markersize=6)
    plt.scatter(pressures, standardDevs, label='500 & 50 V/cm')
    plt.axvline(x=745, color='grey', linestyle='--')
    plt.text(1250, 3, 'Right = 50 V/cm Left= 500 V/cm', bbox=dict(facecolor='white', edgecolor='white', boxstyle='round'))
    plt.xlabel('Pressure (Torr)')
    plt.ylabel('$\sigma$ (mm/$\sqrt{\mathrm{mm}}$)')
    plt.ylim(bottom=0)
    if meanFrom == 'DataSet' :
        plt.title('Sigma From Mean Calculated Per dataset')
    else:
        plt.title('Sigma From Mean ' + str(meanFrom))
    plt.savefig('SigmaData')
    plt.close()

def calibrationProcessing():
    filepaths_list = []

    # Get the current directory where this script is located
    current_dir = os.path.dirname(os.path.abspath(__file__))
    dirList=[]

    dirList.append(os.path.join(current_dir,'Calibration11142023'))
    dirList.append(os.path.join(current_dir,'Calibration11152023'))
    dirList.append(os.path.join(current_dir,'Calibration11162023'))

    chargePerChPerCalibration=[]
    for dir in dirList:
        filepaths_list = []

        filepaths_list = [os.path.join(dir, file) for file in os.listdir(dir) if file.endswith(".txt") and "nA" in file]
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


        fig, axs = plt.subplots(3, 5, figsize=(15, 9))

        # Flatten the axs array for easy iteration
        axs = axs.flatten()

        chargePerReset=[] #14 charges for each channel for this directory
        for i in range(14):
            y = [0.01, 0.07, 0.1, 0.5]
            x = [1 / averageRtdLists[0][i], 1 / averageRtdLists[1][i], 1 / averageRtdLists[2][i], 1 / averageRtdLists[3][i]]

            # Perform linear regression
            coefficients = np.polyfit(x, y, 1) #returns slope and y intercept of fit as coefficients[0] and coefficients[1]
            polynomial = np.poly1d(coefficients)
            line_of_best_fit = polynomial(x)

            chargePerReset.append(coefficients[0])
            # Plot the data points
            axs[i].plot(x, y, 'o', label=f'ChargePerReset: {coefficients[0]:.4f}')
            axs[i].plot(x, y, 'o', label='Data points')

            # Plot the line of best fit
            axs[i].plot(x, line_of_best_fit, label=f'Line of best fit: {coefficients[0]:.4f}x + {coefficients[1]:.4f} ')


            # Set labels and title for the subplot
            axs[i].set_xlabel('1/(Avg Rtd) (1/s)')
            axs[i].set_ylabel('Input I (A)')
            axs[i].set_title('Ch ' + str(i) + '  1/(RTD) per Input I')

            # Show legend with smaller font size
            axs[i].legend(fontsize=6)

        # Adjust layout for better spacing
        plt.tight_layout()

        fig.tight_layout()
        fig.savefig(str(dir.split('/')[-1])+'.png', dpi=500)

        chargePerChPerCalibration.append(chargePerReset)
    return chargePerChPerCalibration

def findRingMeans(chstart = 0):
    # Wellseley areas, the start and stop of each channel where edges are between concentric rings
    minmax = [[0.000, 0.670], [0.670, 1.386], [1.386, 2.106], [2.106, 2.981],
              [2.981, 4.005], [4.005, 4.994], [4.994, 6.015], [6.015, 7.481],
              [7.481, 9.994], [9.994, 12.497], [12.497, 15.023], [15.023, 19.996],
              [19.996, 24.962], [24.962, 30.026], [30.026, 39.977], [39.977, 50.065]]


    ring_avgs = [round((2/3)*(((end**3)-(start**3))/((end**2)-(start**2))), 3) for start, end in minmax]

    # Fiducialize
    ring_avgs = ring_avgs[chstart:]

    return ring_avgs

def normalize_data(data, total_max):
    """Normalizes to a given maxiumum value, denoted by the tallest diffusion peak"""
    max_value = max(data)
    if max_value != 0:
        normalized_data = [value * (total_max / max_value) for value in data]
    else:
        # Handle the case when max_value is zero
        normalized_data = [0.0] * len(data)
    return normalized_data


def calc_quartiles_bounds(data):
    """Box and Whisker Analysis"""
    # Calculate the 25th percentile (Q1) of the 'data'
    q1 = np.percentile(data, 25)

    # Calculate the median (50th percentile or Q2) of the 'data'
    q2 = np.percentile(data, 50)

    # Calculate the 75th percentile (Q3) of the 'data'
    q3 = np.percentile(data, 75)

    # Calculate the interquartile range (IQR) by subtracting Q1 from Q3
    iqr = q3 - q1

    # Calculate the upper bound for outliers using 1.5 times the IQR above Q3
    upper_bound = q3 + 1.5 * iqr

    # Calculate the lower bound for outliers using 1.5 times the IQR below Q1
    lower_bound = q1 - 1.5 * iqr

    return q1, q2, q3, upper_bound, lower_bound


def weighted_average_rtd_bins(rtd_list, num_bins=50):
    """Calculates a weighted average based on the number of rtds in each bin, using a default 50 bins"""
    if len(rtd_list) == 0:
        weighted_mean = 0
    else:
        # Calculate the bin edges
        bin_edges = np.linspace(min(rtd_list), max(rtd_list), num_bins + 1)

        # Bin the values and get the bin counts
        binned_values, _ = np.histogram(rtd_list, bins=bin_edges)

        # Get the bin indices for each value in rtd_list
        bin_indices = np.digitize(rtd_list, bin_edges) - 1

        # Multiply each value by its weight
        weighted_values = [rtd_list[i] * binned_values[bin_indices[i]] if bin_indices[i] < num_bins else 0 for i in
                           range(len(rtd_list))]

        # Calculate the weighted sum of 'weighted_values'
        weighted_sum = sum(weighted_values)

        # Calculate the sum of weights, which involves summing values from 'binned_values' based on 'bin_indices'
        # If 'bin_indices' is out of bounds (greater than or equal to 'num_bins'), use 0 for that index
        sum_weights = sum(
            [binned_values[bin_indices[i]] if bin_indices[i] < num_bins else 0 for i in range(len(rtd_list))])

        # Calculate the weighted mean by dividing the weighted sum by the sum of weights
        weighted_mean = weighted_sum / sum_weights

    return weighted_mean

def zscoreAnalysis(data):
    mean = np.mean(data)
    std = np.std(data)
    data_zscore = [abs((point - mean)/std) for point in data]
    return data_zscore

def areaNormalize(totalQ_PerChannel,  chstart=0):
    """Normalizes the y axis to area"""

    # Calculate the areas and get the errors, but pass the starting channel in case you have fiducialized
    areas, areaErr = calculateAreas(chstart)

    # Area Normalize
    Q_perArea = [totalQ / a for totalQ, a in zip(totalQ_PerChannel, areas)]

    return Q_perArea, areaErr


def calculateAreas(chstart = 0):
    """Calculates the areas per copper ring, with the start/stop point being between consecutive rings"""

    # Wellseley areas, the start and stop of each channel where edges are between concentric rings
    minmax = [[0.000, 0.670], [0.670, 1.386], [1.386, 2.106], [2.106, 2.981],
              [2.981, 4.005], [4.005, 4.994], [4.994, 6.015], [6.015, 7.481],
              [7.481, 9.994], [9.994, 12.497], [12.497, 15.023], [15.023, 19.996],
              [19.996, 24.962], [24.962, 30.026], [30.026, 39.977], [39.977, 50.065]]

    # Area based on outer radius - Area based on inner radius for each channel
    ring_areas = [round((3.14159 * (end ** 2)) - (3.14159 * (start ** 2)), 3) for start, end in minmax]

    #Fiducialize
    ring_areas = ring_areas[chstart:]

    # Error on Q/A is in the radius measurement Q/(pi*R1^2 - pi*R2^2)
    # Error Propagation:
    # delta(Q/A) = sqrt([d/dR1{Q/pi*R1^2}]^2 * sigmaR1^2   +    [d/dR2{-Q/pi*R2^2}]^2 * sigmaR2^2)
    # delta(Q/A) = sqrt([-2Q/pi*R1^3]^2 * sigmaR1^2   +    [2Q/pi*R2^3]^2 * sigmaR2^2)
    # delta(Q/A) = (2Q/pi) * sqrt((sigmaR1^2*R2^6 + sigmaR2^2*R1^6) / (R1^6 * R2^6))
    # delta(Q/A) = (2Q*sigmaR/pi*R1^3*R2^3) * sqrt(R2^6 + R1^6))

    # Sigma radius is error on gerber file ruler (statistical) 0.001 mm
    # We dont have Q here, so just do (2*sigmaR/pi*R1^3*R2^3) * sqrt(R2^6 + R1^6)) and multiply by Q later in plot def
    deltaAErr = [(2*0.001*np.sqrt((start**6)+(end**6)))/(3.14159*(start**3)*(end**3)) for start, end in minmax]

    return ring_areas, deltaAErr


def plot(data_dict,  errorsPerPressureAbove,errorsPerPressureBelow,  save_file_name, title, fiducialize_num_start=0, fiducialize_num_end=16,
         area_normalize=True, plot_charge=True):
    """Input a dictionary with pressure/field keys"""

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'gold', 'lime', 'pink', 'silver', 'gray', 'maroon', 'olive', 'teal',
              'navy']

    # Define a base list of channel numbers
    ch = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

    # Define a base list of mm measurments, where the value here is the midpoint radius of each channel
    #mm = [0.335, 1.028, 1.746, 2.5435, 3.493, 4.4995, 5.5045, 6.748, 8.7375, 11.2455, 13.76, 17.5095, 22.479, 27.494,34.0015, 44.021]
    mm = findRingMeans()
    # Fiducialize the channel and mm
    mm = mm[:fiducialize_num_end]
    mm = mm[fiducialize_num_start:]
    ch = ch[:fiducialize_num_end]
    ch = ch[fiducialize_num_start:]

    fig, ax = plt.subplots()
    n = 0

    chargeList=[]
    errsAboveList=[]
    errsBelowList=[]

    for key, value in data_dict.items():

        # Get the pressure value from the key
        pressure=key

        plotting_data = value

        # Fiducialize the errors, remember we have a list of 16 channels for each pressure so pick the current list [n]
        yerrAbove = errorsPerPressureAbove[n][:fiducialize_num_end]
        yerrAbove = yerrAbove[fiducialize_num_start:]
        yerrBelow = errorsPerPressureBelow[n][:fiducialize_num_end]
        yerrBelow = yerrBelow[fiducialize_num_start:]

        if area_normalize:

            # Area normalize everything, and get the errors on the area while doing so
            plotting_data, areaErr = areaNormalize(value, fiducialize_num_start)

            # Area normalize your errors as well
            yerrAbove, areaErr = areaNormalize(yerrAbove)
            yerrBelow, areaErr = areaNormalize(yerrBelow)

            # The area error is always the same, but you need to do Q/areaError
            # areaErr has been set up for this already by doing 1/err, now just multiply by charge
            areaErrFinal=[charge * areaer for charge, areaer in zip(plotting_data,areaErr)]

            # Add the area error above and below each data point
            yerrAbove = [yerrAb + areaEr for yerrAb,areaEr in zip(yerrAbove,areaErrFinal)]
            yerrBelow = [yerrbel + areaEr for yerrbel,areaEr in zip(yerrBelow,areaErrFinal)]

            if plot_charge:
                ax.set_ylabel('Total Charge (C) per mm\u00b2')
            if not plot_charge:
                ax.set_ylabel('Number of RTDs per mm^2')

        if not area_normalize:

            plotting_data = value

            if plot_charge:
                ax.set_ylabel('Total Charge (C)')
            if not plot_charge:
                ax.set_ylabel('Number of RTDs')


        #chargeList.append('Pressure '+ str(round(float(pressure) * 51.715, 1)) + ': ')
        #errsAboveList.append('Pressure '+ str(round(float(pressure) * 51.715, 1)) + ': ')
        #errsBelowList.append('Pressure '+ str(round(float(pressure) * 51.715, 1)) + ': ')

        #chargeList.append(plotting_data)
        #errsAboveList.append(yerrAbove)
        #errsBelowList.append(yerrBelow)
        pressure = round(float(pressure) * 51.715, 1)
        formatted_pressure = f'{pressure:06.1f}'

        ax.plot(mm, plotting_data, marker='.', label=str(formatted_pressure) + " Torr", color=colors[n])
        ax.errorbar(mm, plotting_data, yerr=(yerrBelow,yerrAbove), fmt='none', capsize=6, markersize='6', color=colors[n])
        n += 1

    handles, labels = ax.get_legend_handles_labels()
    # sort both labels and handles by labels
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    ax.legend(handles, labels)

    ax.set_xlabel('Radius (mm)')
    ax.set_title(title)

    fig.savefig(save_file_name)


def get_filepaths(current_directory=True):
    """input a text file containing a 2d list with 16 lists of rtds, parsed from
    the root file and converted to seconds. Specify if your data directory is in this directory or not"""
    filepaths_list = []
    if current_directory:
        # Get the current directory where this script is located
        current_dir = os.path.dirname(os.path.abspath(__file__))
        calibratedPressureFileNames = []
        dirList = []

        dirList.append(os.path.join(current_dir, 'Calibration11142023'))
        dirList.append(os.path.join(current_dir, 'Calibration11152023'))
        dirList.append(os.path.join(current_dir, 'Calibration11162023'))

        for dir in dirList:
            filepaths_list = [os.path.join(dir, file) for file in os.listdir(dir) if
                              file.endswith(".txt") and "kV" in file]
            for file in filepaths_list:
                calibratedPressureFileNames.append(file)

    else:
        print("update get_filepaths definition with code to get to your directory")

    return calibratedPressureFileNames


def make_data_dictionary(chargePerResetPerCh, filepaths, fiducialize_num_start=0, fiducialize_num_end=16, background_subtracted=True,
                         charge_convert=True):
    """For the charge conversion, if it is true, make sure the vdd values are saved in the same directory
     as the run1_times.txt file, and in an excel sheet titled 'VddBeforeAfter.xlsx' with the Vdd values in sequential order
      from 1-16 in the second column"""
    errsAbove = []
    errsBelow = []
    data_dictionary = {}
    background_dictionary = {}
    num = 0


    # Go through each individual pressure run and do analysis
    for filepath in filepaths:
        num += 1

        # Initialize a list to hold the 16 lists of reset time differences for data
        reset_time_diffs = [[] for _ in range(16)]

        # Initialize a list to hold the 16 lists of reset time differences for backgrounds
        background_reset_time_diffs = [[] for _ in range(16)]
        background_data = None

        # Read in the data from the text file and populate the list
        with open(filepath, "r") as f:
            # Read the contents of the file and use 'ast.literal_eval()' to safely parse it as a Python literal
            data = ast.literal_eval(f.read())

            # Iterate over the elements of 'data' using enumeration
            maxTimeStamp=0
            for i, timeStamps in enumerate(data):

                # Check if there are rtds, if not leave the list empty as it already is
                if len(timeStamps) > 0:
                    # Filter out values in 'rtd_list' that are less than 0.05 because they are considered not physical
                    # rtd=0.05s means I=CV*4/rtd=(1.0E-11)*0.4v/0.05=8e-11A which is 80 pico amps, and we saw max 60 picoamps on the single channel
                    # we saw huge spikes below this threshold and thought it was due to something in the DAQ
                    rtd_list = []
                    timeStamps= [timestamp * 200 / 30e6 for timestamp in timeStamps]
                    if timeStamps[-1] > maxTimeStamp:
                        maxTimeStamp = timeStamps[-1]
                    for index in range(0,len(timeStamps)):
                        if index !=0:
                            rtd = timeStamps[index]-timeStamps[index - 1]
                            rtd_list.append(rtd)
                    #varghese already cuts in his analysis
                    #reset_time_diffs[i] = [value for value in rtd_list if value >= 0.0088]
                    reset_time_diffs[i] = [value for value in rtd_list]

                else:
                    continue

        # Do quartile analysis to remove outliers (not doing this really hurts the error bars)
        #fig, axs = plt.subplots(4, 4, figsize=(12, 12))
        #axs = axs.ravel()
        meanRtdPerCh=[0]*(len(reset_time_diffs))

        # Parse the pressure and feild from the headers
        if "Hgbelowatm" in filepath:
            inHg = filepath.split('inHgbelowatm')[0].split('100Hz_')[1]
            inHg = 29.9213 - float(inHg)  # below 1 atm (1atm=29.9213 inHg)
            pressure = inHg / 2.036  # (in psi)

        elif "psi" in filepath:
            RelativePsi = filepath.split('psi')[0].split('100Hz_')[1]
            pressure = 14.6959 + float(RelativePsi)
        else:
            print("error in parsing file headers")
            print(filepath)

        fig, axs = plt.subplots(4, 4, figsize=(12, 12))
        axs = axs.ravel()
        for i in range(0, len(reset_time_diffs)):

            # Need to make this greater than one reset because otherwise there is a division by 0 error
            if len(reset_time_diffs[i]) > 1:
                #q1, q2, q3, upper_bound, lower_bound = calc_quartiles_bounds(reset_time_diffs[i])

                #get z scorea
                data_zscores = zscoreAnalysis(reset_time_diffs[i])

                #filter based on 2 sigma
                reset_time_diffs[i] = [x for x, y in zip(reset_time_diffs[i], data_zscores) if (y < 3)]
                #filteredzscores= [y for x, y in zip(reset_time_diffs[i], data_zscores) if (y < 2)]

                #iterate twice
                data_zscores = zscoreAnalysis(reset_time_diffs[i])
                reset_time_diffs[i] = [x for x, y in zip(reset_time_diffs[i], data_zscores) if (y < 2)]

                # Get the new mean
                meanRtd = np.mean(reset_time_diffs[i])
                meanRtdPerCh[i] = meanRtd

                # Uncomment to take a look at your rtds for each individual channel on a different plot per pressure

                '''axs[i].hist(reset_time_diffs[i], bins=50) # loook into this function for bins
                axs[i].axvline(meanRtd, color='blue', linestyle='--')
                axs[i].set_xlabel("rtd (s)")
                axs[i].set_ylabel("Frequency")
                axs[i].legend(fontsize='6')
                axs[i].set_title("Ch{}".format(i + 1))'''

        #fig.savefig(str(round(float(pressure) * 51.715, 1))+'.png')

        # Output a list of total charge seen per channel
        if charge_convert:
            if "11142023" in filepath: #figure out which calibration charge value to use
                num=0
            if "11152023" in filepath:
                num=1
            if "11162023" in filepath:
                num=2
            deltaQ_perReset_perChannel = chargePerResetPerCh[num] #grab that list of 14 charge values

            totalQ_PerChannel = []
            errorsPerPressureAbove = []
            errorsPerPressureBelow = []

            # Only do this for the fiducial volume
            for i in range(fiducialize_num_start, fiducialize_num_end):

                #number of resets * deltaChargePerReset = total charge seen on that channel
                if np.average(reset_time_diffs[i]) ==0:
                    totalQ_PerChannel.append(0)
                else:
                    totalQ_PerChannel.append(deltaQ_perReset_perChannel[i]*maxTimeStamp/np.average(reset_time_diffs[i]) )

                # -------------------- Make your errors here ----------------------
                if len(reset_time_diffs[i]) == 0:
                    errorsPerPressureAbove.append(0)
                    errorsPerPressureBelow.append(0)
                else:
                    # Determine your desired confidence interval
                    ci = 1

                    # Find your error on total Q seen:
                    # I = Q/t = CV*4/rtd
                    # Q = qt/rtd
                    # Now propagate the errors on rtds:
                    # SigmaQ = sqrt[(d/dRTD[qt/meanRTD])^2 * rtdSigma^2]
                    # SigmaQ = (qt/meanRTD^2) * rtdSigma
                    systematics = ci * np.std(reset_time_diffs[i]) * deltaQ_perReset_perChannel[i] * 1800/((meanRtdPerCh[i])**2)

                    # Add the the one quanta error, but only to the top (You can miss extra reset if you stop in middle)
                    # Add another quanta for when you make the rtd list if there is an odd number
                    errorsPerPressureAbove.append(systematics + 2* deltaQ_perReset_perChannel[i])
                    errorsPerPressureBelow.append(systematics)

        # Otherwise, just do the total number of rtds
        if not charge_convert:
            totalQ_PerChannel = [len(rtd_list) for rtd_list in reset_time_diffs]

        # For each pressure add your 16 channel list of errors
        errsAbove.append(errorsPerPressureAbove)
        errsBelow.append(errorsPerPressureBelow)


        # Add the data for this pressure/field to your dictionary with a key

        data_dictionary[pressure] = totalQ_PerChannel

    return data_dictionary,  errsAbove,errsBelow

chargePerResetPerCh= calibrationProcessing()
filepaths = get_filepaths()
data_dict, errorsPerPressureAbove,errorsPerPressureBelow = make_data_dictionary(chargePerResetPerCh,filepaths, fiducialize_num_start=0,
                                                                     fiducialize_num_end=11, background_subtracted=False,
                                                                     charge_convert=True)

plot(data_dict, errorsPerPressureAbove,errorsPerPressureBelow, 'ZscoreRTDOutlierMethod5.png', 'Charge Diffusion (Zscore> 3 on iteration 1, Zscore> 2 on iteration 2)',
     fiducialize_num_start=0, fiducialize_num_end=11, area_normalize=True, plot_charge=True)

# fixed mean is based on spot orientation analysis posted in slack
plotDt(data_dict,  errorsPerPressureAbove,errorsPerPressureBelow ,fiducialize_num_start=0, fiducialize_num_end=11, DtUnits='cm^s/s', area_normalize=True, mean_from_data=False, fixed_mean = 0.01 )
