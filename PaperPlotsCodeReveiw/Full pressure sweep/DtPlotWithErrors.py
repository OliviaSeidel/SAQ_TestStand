import matplotlib.pyplot as plt
import ast
import numpy as np
import os
import pandas as pd

def sum_list_values(lst):
    total = 0
    for value in lst:
        total += value
    return total

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

    # Area based on outter radius - Area based on inner radius for each channel
    ring_areas = [round((3.14159 * (end ** 2)) - (3.14159 * (start ** 2)), 3) for start, end in minmax]

    #Fiducialize
    ring_areas = ring_areas[chstart:]

    # x/area = c/mm^2
    # error on this is charge/a+a(0.001) - charge/a-a(0.001) = charge *(1/a+a(0.001) - 1/a-a(0.001))
    # sigma area error is just sigma-radius converted to area. Sigma radius is error on gerber file ruler (statistical)
    deltaA=3.14159 * ((0.001) ** 2)

    # Error in area normalization:
    deltaAErr = [(1/(r_area+deltaA))-(1/(r_area-deltaA)) for r_area in ring_areas]

    return ring_areas, deltaAErr


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

    # Drift velocities pulled from Varghese simulation plots
    Vd = [52, 50, 48, 45, 43, 39, 35, 37, 40, 43, 45, 46, 47, 48, 49]

    for key, value in data_dict.items():
        plotting_data = None

        pressure = key.split('psi')[0].replace('neg', '-').replace('.os', '')
        pressure = round(float(pressure) * 51.715, 1)
        pressures.append(pressure)

        if area_normalize:
            plotting_data, area_err = areaNormalize(value, fiducialize_num_start)

        if not area_normalize:
            plotting_data = value

        plotting_sets.append(plotting_data)
        i += 1

    if mean_from_data:
        totChargePerPress, standardDevs, mean, numerator = calculate_standard_deviations_mean(mm, plotting_sets, meanVal = 'DataSet', DtUnits=DtUnits)
        plot_sigma(pressures, standardDevs, numerator, totChargePerPress, errsAbove, errsBelow, DtUnits = DtUnits,meanFrom= 'DataSet')
        plot_diffusion_coefficient(pressures, numerator, errsAbove, errsBelow, Vd, standardDevs, 'DtData', meanFrom='DataSet', DtUnits= DtUnits)
    if not mean_from_data:
        mean = fixed_mean
        totChargePerPress, standardDevs, mean, numerator = calculate_standard_deviations_mean(mm, plotting_sets, meanVal = mean, DtUnits=DtUnits)
        plot_sigma(pressures, standardDevs,  numerator, totChargePerPress, errsAbove, errsBelow,DtUnits = DtUnits,meanFrom = mean)
        plot_diffusion_coefficient(pressures, numerator, totChargePerPress,errsAbove, errsBelow, Vd, standardDevs, 'DtData',meanFrom = mean, DtUnits = DtUnits)


    return standardDevs

def get_fitting_data(mean, mm, plotting_data):
    if mean > mm[1]:
        fittingmm = mm[2:]
        fittingdata = plotting_data[2:]
    elif mean < mm[1] and mean > mm[0]:
        fittingmm = mm[1:]
        fittingdata = plotting_data[1:]
    else:
        fittingmm = mm
        fittingdata = plotting_data
    return fittingmm, fittingdata


def calculate_standard_deviations_mean(mm, plotting_sets, meanVal= 'DataSet', DtUnits='cm^s/s'):
    standardDevs = []
    for i in range(0, len(plotting_sets)):
        ydata_norm = np.array(plotting_sets[i])
        bins = 300

        counts, bin_edges = np.histogram(mm, bins=bins, weights=ydata_norm)
        non_zero_bins = counts > 0
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_centers = bin_centers[non_zero_bins]
        #counts = list(filter(lambda x: x != 0, counts))
        # Calculate the weighted sum of bin_centers multiplied by counts
        weighted_sum = np.sum(bin_centers * counts[non_zero_bins])

        # Calculate the total sum of counts
        total_count = np.sum(counts)

        # Calculate the mean
        if meanVal == 'DataSet':
            mean = weighted_sum / total_count
            print("Mean of histogrammed data: ", mean)
        else:
            mean = meanVal
            print("Set Mean: ", mean)

        centerminusmean = ((bin_centers - mean) ** 2)
        countstimesdiff = counts[non_zero_bins] * centerminusmean
        totChargePerPress=sum_list_values(counts[non_zero_bins])
        nMinus1 = (totChargePerPress)
        numerator = sum_list_values(countstimesdiff)
        std_mm = np.sqrt(numerator / nMinus1)
        standardDevs.append(std_mm)

    return totChargePerPress, standardDevs,mean, numerator

def plot_sigma(pressures, standardDevs,  numerator, totChargePerPress, errsAbove, errsBelow, DtUnits = 'DtUnits',meanFrom= 'DataSet'):
    YerrsAbove = []
    YerrsBelow = []

    for i, std_mm in enumerate(standardDevs):

        errAbove = (numerator * np.sum(errsAbove[i])) / (20 * np.sqrt(numerator/totChargePerPress) * (totChargePerPress ** 2))
        errBelow = (numerator * np.sum(errsBelow[i])) / (20 * np.sqrt(numerator/totChargePerPress) * (totChargePerPress ** 2))

        YerrsAbove.append(errAbove)
        YerrsBelow.append(errBelow)
    standardDevs = [x/np.sqrt(100) for x in standardDevs] #100mm is 10 cm for 10 cm drift
    plt.figure()
    plt.errorbar(pressures, standardDevs, yerr=(YerrsBelow, YerrsAbove), fmt='none', capsize=6, markersize=6)
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


def plot_diffusion_coefficient(pressures, numerator, totChargePerPress, errsAbove, errsBelow, Vd, standardDevs, filename, meanFrom = 'DataSet', DtUnits= 'DtUnits'):
    plt.figure()
    Dtvals = []
    YerrsAbove = []
    YerrsBelow = []

    for i, std_mm in enumerate(standardDevs):
        Dt = (((std_mm * 0.1) ** 2) * Vd[i] * 0.1) / (2 * 11 * (1e-6))
        Dtvals.append(Dt)

        #errAbove =  2 * np.sqrt(totChargePerPress/ np.sum(errsAbove[i]))
        #errBelow =  2 * np.sqrt(totChargePerPress / np.sum(errsBelow[i]))
        errAbove = (numerator * ((0.1)**3) * Vd[i] * np.sum(errsAbove[i])) / ((22e-6) * (totChargePerPress**2))
        errBelow =(numerator * ((0.1)**3) * Vd[i] * np.sum(errsBelow[i])) / ((22e-6) * (totChargePerPress**2))

        YerrsAbove.append(errAbove)
        YerrsBelow.append(errBelow)

    plt.scatter(pressures, Dtvals, label='500 & 50 V/cm')
    plt.text(1250, 17500, 'Right = 50 V/cm Left= 500 V/cm',
             bbox=dict(facecolor='white', edgecolor='white', boxstyle='round'))
    plt.axvline(x=745, color='grey', linestyle='--')
    plt.errorbar(pressures, Dtvals, yerr=(YerrsBelow, YerrsAbove), fmt='none', capsize=6, markersize=6)
    plt.xlabel('Pressure (Torr)')
    plt.ylabel('Dt (cm^2/s)')
    plt.ylim(bottom=0)
    if meanFrom == 'DataSet' :
        plt.title('Dt From Mean Calculated Per dataset at 500  V/cm')
    else:
        plt.title('Dt From Mean ' + str(meanFrom) + ' at 500 V/cm')

    plt.savefig(filename)
    plt.close()

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

def plot(data_dict,  errorsPerPressureAbove,errorsPerPressureBelow,  save_file_name, title, fiducialize_num_start=0, fiducialize_num_end=16,
         area_normalize=True, plot_charge=True):
    """Input a dictionary with pressure/field keys"""

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'gold', 'lime', 'pink', 'silver', 'gray', 'maroon', 'olive', 'teal',
              'navy']

    # Define a base list of channel numbers
    ch = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

    # Define a base list of mm measurments, where the endpoint of each channel is between two concentric rings
    mm = [0.335, 1.028, 1.746, 2.5435, 3.493, 4.4995, 5.5045, 6.748, 8.7375, 11.2455, 13.76, 17.5095, 22.479, 27.494,
          34.0015, 44.021]

    # Fiducialize the channel and mm
    mm = mm[:fiducialize_num_end]
    mm = mm[fiducialize_num_start:]
    ch = ch[:fiducialize_num_end]
    ch = ch[fiducialize_num_start:]

    totYerrAbove=[]
    totYerrBelow = []


    fig, ax = plt.subplots()
    n = 0
    for key, value in data_dict.items():

        # Get the pressure value from the key
        pressure = key.split('psi')[0].replace('neg', '-').replace('.os', '')

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
            areaErrFinal = [charge * areaer for charge, areaer in zip(plotting_data, areaErr)]

            # Add the area error above and below each data point
            yerrAbove = [yerrAb + areaEr for yerrAb, areaEr in zip(yerrAbove, areaErrFinal)]
            yerrBelow = [yerrbel + areaEr for yerrbel, areaEr in zip(yerrBelow, areaErrFinal)]

            # If the errror dips to less than 0 charge it doesnt make sense, so set it to 0
            if any(value < 0 for value in yerrBelow) or any(value < 0 for value in yerrAbove):
                yerrBelow = [value if value >= 0 else 0 for value in yerrBelow]
                yerrAbove = [value if value >= 0 else 0 for value in yerrAbove]

            totYerrAbove.append(yerrAbove)
            totYerrBelow.append(yerrBelow)

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

        # Apply your filters if you only want certain pressures:
        if round(float(pressure) * 51.715, 1) > 720:

            ax.plot(mm, plotting_data, marker='.', label=str(round(float(pressure) * 51.715, 1)) + ' Torr', color=colors[n])
            ax.errorbar(mm, plotting_data, yerr=(yerrBelow,yerrAbove), fmt='none', capsize=6, markersize='6', color=colors[n])
        n += 1

    ax.set_xlabel('Radius (mm)')
    ax.set_title(title)
    ax.legend()
    fig.savefig(save_file_name)

    # return the total charge errors above and below for dt analysis
    return totYerrAbove, totYerrBelow



def get_filepaths(textfile, current_directory=True):
    """input a text file containing a 2d list with 16 lists of rtds, parsed from
    the root file and converted to seconds. Specify if your data directory is in this directory or not"""
    filepaths_list = []
    if current_directory:
        # Get the current directory where this script is located
        current_dir = os.path.dirname(os.path.abspath(__file__))

        filepaths_list = []

        # Recursively traverse through the directory tree starting from current_dir
        for root, dirs, files in os.walk(current_dir):
            # Iterate through the subdirectories in the current directory
            for dir_name in dirs:
                # Create the full path to the subdirectory
                dir_path = os.path.join(root, dir_name)

                # Create a full path by appending 'textfile' to the subdirectory path
                dir_path = os.path.join(dir_path, textfile)

                # Append the full path to the filepaths_list
                filepaths_list.append(dir_path)

        # Sort the list of file paths alphabetically
        filepaths_list = sorted(filepaths_list)
    else:
        print("update get_filepaths definition with code to get to your directory")

    return filepaths_list


def make_data_dictionary(filepaths, fiducialize_num_start=0, fiducialize_num_end=16, background_subtracted=True,
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
            for i, rtd_list in enumerate(data):
                # Check if there are rtds, if not leave the list empty as it already is
                if len(rtd_list) > 0:
                    # Filter out values in 'rtd_list' that are less than 0.05 because they are considered not physical
                    # rtd=0.05s means I=CV*4/rtd=(1.0E-11)*0.4v/0.05=8e-11A which is 80 pico amps, and we saw max 60 picoamps on the single channel
                    # we saw huge spikes below this threshold and thought it was due to something in the DAQ
                    reset_time_diffs[i] = [value for value in rtd_list if value >= 0.05]

                else:
                    continue

        # Do quartile analysis to remove outliers (not doing this really hurts the error bars)
        fig, axs = plt.subplots(4, 4, figsize=(12, 12))
        axs = axs.ravel()
        meanRtdPerCh=[0]*(len(reset_time_diffs))

        for i in range(0, len(reset_time_diffs)):

            # Need to make this greater than one reset because otherwise there is a division by 0 error
            if len(reset_time_diffs[i]) > 1:
                q1, q2, q3, upper_bound, lower_bound = calc_quartiles_bounds(reset_time_diffs[i])

                # Discard the outliers beyond the upper and lower bound (see "Test Stand Update May 26" for more info)
                reset_time_diffs[i] = [x for x in reset_time_diffs[i] if (x <= (upper_bound) and x >= lower_bound)]

                # Get the weighted mean RTD
                meanRtd=weighted_average_rtd_bins(reset_time_diffs[i])
                meanRtdPerCh[i] = meanRtd
                # Uncomment to take a look at your rtds for each individual channel on a different plot per pressure
                '''
                axs[i].hist(reset_time_diffs[i], bins=50)
                axs[i].axvline(q1, color='blue', linestyle='--', label='q1')
                axs[i].axvline(q2, color='purple', linestyle='--', label='q2')
                axs[i].axvline(q3, color='green', linestyle='--', label='q3')
                axs[i].axvline(upper_bound, color='black', linestyle='--', label='upper bound')
                axs[i].axvline(lower_bound, color='yellow', linestyle='--', label='lower bound')
                # axs[i].axvline(avgrtd, color='blue', linestyle='--')
                axs[i].set_xlabel("rtd")
                axs[i].set_ylabel("Frequency")
                axs[i].legend(fontsize='6')
                axs[i].set_title("Ch{}".format(i + 1))
                '''
        # fig.savefig(str(num))

        # Get your vdd file names
        filepathVdds = filepath.replace('run1_times.txt', 'VddBeforeAfter.xlsx')

        # Read the Excel file into a DataFrames
        df = pd.read_excel(filepathVdds)

        # Extract values from the second column and convert them to a list
        vdd_list = df.iloc[:, 1].tolist()


        # Get background data if needed
        if background_subtracted:
            filepathBackground = filepath.replace('run1', 'background')
            with open(filepathBackground, "r") as f:
                data = ast.literal_eval(f.read())
                for i, rtd_list_background in enumerate(data):
                    if len(rtd_list_background) > 0:
                        background_reset_time_diffs[i] = [value for value in rtd_list_background if
                                                          value >= 0.05]  # throw out values less than this because not physical
            # output a list of total charge seen per channel
            totalQ_PerChannel_background = []

            for i in range(0, 16):
                totalQ_PerChannel_background.append(len(background_reset_time_diffs[i]) * deltaQ_perReset_perChannel[i])

            # Feducial volume cut:
            totalQ_PerChannel_background = totalQ_PerChannel_background[:fiducialize_num_end]
            background_data = totalQ_PerChannel_background[fiducialize_num_start:]

        # Output a list of total charge seen per channel
        if charge_convert:

            # Convert the total rtds to total charge by doing cv*num_rtds which is also deltaQ*numResets=total charge
            # Assume c=10pF (we measured this and it was very close on all channels), take vdd and multiply times 4, convert to V
            deltaQ_perReset_perChannel = [1.0e-11 * vdd * 4 * 1e-3 for vdd in vdd_list]

            totalQ_PerChannel = []
            errorsPerPressureAbove = []
            errorsPerPressureBelow = []

            # Only do this for the fiducial volume
            for i in range(fiducialize_num_start, fiducialize_num_end):

                #number of resets * deltaChargePerReset = total charge seen on that channel
                totalQ_PerChannel.append(len(reset_time_diffs[i]) * deltaQ_perReset_perChannel[i])

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
                    # SigmaQ = (qt/meanRTD^2) * rtdSigma
                    systematics = ci * np.std(reset_time_diffs[i]) * deltaQ_perReset_perChannel[i] * 1800/((meanRtdPerCh[i])**2)

                    # Add the the one quanta error, but only to the top (You can miss extra reset if you stop in middle)
                    errorsPerPressureAbove.append(systematics + deltaQ_perReset_perChannel[i])
                    errorsPerPressureBelow.append(systematics)

            if background_subtracted:
                # Subtract background
                totalQ_PerChannel = [(totQ - (len(background_rtd) * deltaQ)) for totQ, background_rtd, deltaQ in
                                     zip(totalQ_PerChannel, background_reset_time_diffs, deltaQ_perReset_perChannel)]


        # Otherwise, just do the total number of rtds
        if not charge_convert:
            totalQ_PerChannel = [len(rtd_list) for rtd_list in reset_time_diffs]

        # For each pressure add your 16 channel list of errors
        errsAbove.append(errorsPerPressureAbove)
        errsBelow.append(errorsPerPressureBelow)

        # Parse the pressure and feild from the headers
        pressure = filepath.split('/')[-2].split('_')[-1].replace('psi', '').replace('p', '.')
        field = filepath.split('/')[-2].split('_')[0].replace('vPerCm', 'V/cm')

        # Add the data for this pressure/field to your dictionary with a key
        background_dictionary[pressure + "psi" + field] = background_data
        data_dictionary[pressure + "psi" + field] = totalQ_PerChannel

    return data_dictionary, background_dictionary, errsAbove,errsBelow


filepaths = get_filepaths('run1_times.txt')
data_dict, background_dict, errorsPerPressureAbove,errorsPerPressureBelow = make_data_dictionary(filepaths, fiducialize_num_start=0,
                                                                     fiducialize_num_end=11, background_subtracted=False,
                                                                     charge_convert=True)
errsAbove, errsBelow = plot(data_dict, errorsPerPressureAbove,errorsPerPressureBelow, 'AreaNormalizedChargeWithErrorsAbove1ATM2.png', 'Diffusion at Various Pressures at 500V/cm',
     fiducialize_num_start=0, fiducialize_num_end=11, area_normalize=True, plot_charge=True)

# fixed mean is based on spot orientation analysis posted in slack
plotDt(data_dict, errsAbove, errsBelow,fiducialize_num_start=0, fiducialize_num_end=11, DtUnits='cm^s/s', area_normalize=True, mean_from_data=False, fixed_mean = 0.01 )
