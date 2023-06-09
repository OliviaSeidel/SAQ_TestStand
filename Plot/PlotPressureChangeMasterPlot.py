import matplotlib.pyplot as plt
import ast
import numpy as np
import os
import pandas as pd


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
    q1 = np.percentile(data, 25)
    q2 = np.percentile(data, 50)
    q3 = np.percentile(data, 75)

    iqr = q3 - q1
    upper_bound = q3 + 1.5 * iqr
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

        weighted_sum = sum(weighted_values)
        sum_weights = sum(
            [binned_values[bin_indices[i]] if bin_indices[i] < num_bins else 0 for i in range(len(rtd_list))])
        weighted_mean = weighted_sum / sum_weights

    return weighted_mean


def areaNormalize(totalQ_PerChannel):
    """Normalizes the y axis to area"""
    Q_perArea = [totalQ / a for totalQ, a in zip(totalQ_PerChannel, calculateAreas())]
    return Q_perArea


def calculateAreas():
    """Calculates the areas per copper ring"""
    # Contains the start,end of each copper ring consecutively
    start_and_end_Channel = [0, 0.53, 0.79, 1.3, 1.47, 2.01, 2.19, 2.89, 3.08, 3.93, 4.1, 4.9, 5.1, 5.9, 6.12, 7.37,
                             7.58, 9.92, 10.11, 12.39, 12.59, 14.87, 15.090, 19.8, 20.12, 24.79, 25.12, 29.86, 30.15,
                             39.76, 40.060, 47.7]

    # split into the start and end of each channel
    end_ch, start_ch = start_and_end_Channel[1::2], start_and_end_Channel[::2]

    # make a 2d list of min and max for each channel, where start and end is the middle of the space between each copper ring
    minmax = []
    end = 0  # start at 0mm
    num = 0
    for i in range(0, 16):
        num += 1
        start = end
        if num == 16:
            # If the last channel, count the end as the end
            end = round(end_ch[i], 3)
        else:
            # If not the last channel,count the end as the (end+half the distance to the next ring)
            end = round(end_ch[i] + ((start_ch[i + 1] - end_ch[i]) / 2), 3)  # round to three decimal places
        start_end = [start, end]
        minmax.append(start_end)
    ring_areas = [round((3.14159 * (end ** 2)) - (3.14159 * (start ** 2)), 3) for start, end in minmax]
    return ring_areas


def plot(data_dict, save_file_name, title, fiducialize_num=16, area_normalize=True, peak_normalize=True, plot_charge=True):
    """Input a dictionary with pressure/field keys"""
    # get all of the fields so you can make plots for each field
    field_values = []
    total_max = 0
    for key, value in data_dict.items():
        this_field = key.split('psi')[1].split('V/cm')[0]
        if this_field not in field_values:
            field_values.append(this_field)
        max_value = max(value)
        if max_value > total_max:
            total_max = max_value

    ch = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,12,13,14,15,16]
    # fiducialize after channel specified
    ch = ch[:fiducialize_num]

    # plot on a separate plot the various pressures for each field
    for field_value in field_values:
        fig, ax = plt.subplots()
        for key, value in data_dict.items():
            if field_value in key:
                plotting_data=value
                if area_normalize:
                    if peak_normalize:
                        plotting_data =normalize_data(areaNormalize(value), total_max)
                        if plot_charge:
                            ax.set_ylabel('Charge per mm^2 Peak Normalized')
                        if not plot_charge:
                            ax.set_ylabel('Number of RTDs per mm^2 Peak Normalized')
                    if not peak_normalize:
                        plotting_data = areaNormalize(value)
                        if plot_charge:
                            ax.set_ylabel('Total Charge (C) per mm^2')
                        if not plot_charge:
                            ax.set_ylabel('Number of RTDs per mm^2')
                if not area_normalize:
                    if peak_normalize:
                        plotting_data = normalize_data(value, total_max)
                        if plot_charge:
                            ax.set_ylabel('Charge Peak Normalized')
                        if not plot_charge:
                            ax.set_ylabel('Number of RTDs Peak Normalized')
                    if not peak_normalize:
                        plotting_data = value
                        if plot_charge:
                            ax.set_ylabel('Total Charge (C)')
                        if not plot_charge:
                            ax.set_ylabel('Number of RTDs')
                ax.plot(ch, plotting_data, marker='.', label=pressure)
                ax.set_xlabel('Channel Number')
                ax.set_title(title)
                ax.legend(fontsize=8)
                fig.savefig(save_file_name+field_value)

def get_filepaths(textfile, current_directory=True):
    """input a text file containing a 2d list with 16 lists of rtds, parsed from
    the root file and converted to seconds. Specify if your data directory is in this directory or not"""
    filepaths_list = []
    if current_directory:
        current_dir = os.path.dirname(os.path.abspath(__file__))
    for root, dirs, files in os.walk(current_dir):
        for dir_name in dirs:
            dir_path = os.path.join(root, dir_name)
            dir_path = os.path.join(dir_path, textfile)
            filepaths_list.append(dir_path)
    filepaths_list = sorted(filepaths_list)
    return filepaths_list

def make_data_dictionary(filepaths,  fiducialize_num=16,background_subtracted=True, charge_convert=True):
    """For the charge conversion, if it is true, make sure the vdd values are saved in the same directory
     as the run1.txt file, and in an excel sheet titled 'VddBeforeAfter.xlsx' with the Vdd values in sequential order
      from 1-16 in the second column"""

    data_dictionary = {}
    for filepath in filepaths:
        # Initialize a list to hold the 16 lists of reset time differences
        reset_time_diffs = [[] for _ in range(16)]
        if background_subtracted:
            background_reset_time_diffs = [[] for _ in range(16)]

        # Read in the data from the text file and populate the list
        with open(filepath, "r") as f:
            data = ast.literal_eval(f.read())
            for i, rtd_list in enumerate(data):
                if len(rtd_list) > 0:
                    reset_time_diffs[i] = [value for value in rtd_list if
                                           value >= 0.05]  # throw out values less than this because not physical
        filepathBackground = filepath.replace('run1', 'background')
        with open(filepathBackground, "r") as f:
            data = ast.literal_eval(f.read())
            for i, rtd_list_background in enumerate(data):
                if len(rtd_list_background) > 0:
                    background_reset_time_diffs[i] = [value for value in rtd_list_background if
                                                      value >= 0.05]  # throw out values less than this because not physical
        if charge_convert:
            filepathVdds = filepath.replace('run1.txt', 'VddBeforeAfter.xlsx')
            # Read the Excel file into a DataFrame, skip the first row with the headers
            df = pd.read_excel(filepathVdds, skiprows=1)
            # Extract values from the second column and convert them to a list
            vdd_list = df.iloc[:, 1].tolist()

            # convert the total rtds to total charge by doing cv*num_rtds which is also deltaQ*numResets=total charge
            # assume c=10pF, take vdd and multiply times 4, convert to v
            deltaQ_perReset_perChannel = [1.0e-11 * vdd * 4 * 1e-3 for vdd in vdd_list]

            # output a list of total charge seen per channel
            totalQ_PerChannel = []
            for i in range(0, 16):
                totalQ_PerChannel.append(len(reset_time_diffs[i]) * deltaQ_perReset_perChannel[i])

            # ignore everything past channel 11, this is the fedutial volume cut:
            totalQ_PerChannel = totalQ_PerChannel[:fiducialize_num]
            if background_subtracted:
                # subtract background
                totalQ_PerChannel = [(totQ - (len(background_rtd) * deltaQ)) for totQ, background_rtd, deltaQ in
                                     zip(totalQ_PerChannel, background_reset_time_diffs, deltaQ_perReset_perChannel)]
            data = totalQ_PerChannel
        #Otherwise, just do the total number of rtds
        if not charge_convert:
            data = [len(rtd_list) for rtd_list in reset_time_diffs]

        pressure = filepath.split('/')[-2].split('_')[-1].replace('psi', '').replace('p', '.')
        field = filepath.split('/')[-2].split('_')[0].replace('vPerCm', 'V/cm')

        data_dictionary[pressure + "psi" + field] = data
    return data_dictionary

filepaths = get_filepaths('run1.txt')
data_dict = make_data_dictionary(filepaths, fiducialize_num=11, charge_convert=True)
plot(data_dict,'diffusion','Charge over 5 Minute run', fiducialize_num=11, area_normalize=True, peak_normalize=True, plot_charge=True)
