import matplotlib.pyplot as plt
import ast
import numpy as np
import os
import pandas as pd
from scipy.optimize import curve_fit
from scipy.integrate import quad

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

def areaNormalize(totalQ_PerChannel, chstart=0):
    """Normalizes the y axis to area"""
    Q_perArea = [totalQ / a for totalQ, a in zip(totalQ_PerChannel, calculateAreas(chstart))]
    return Q_perArea

def calculateAreas(chstart = 0):
    """Calculates the areas per copper ring"""
    #wesley Areas
    minmax= [ [ 0.000,  0.670], [ 0.670,  1.386], [ 1.386,  2.106], [ 2.106,  2.981],
           [ 2.981,  4.005], [ 4.005,  4.994], [ 4.994,  6.015], [ 6.015,  7.481],
           [ 7.481,  9.994], [ 9.994, 12.497], [12.497, 15.023], [15.023, 19.996],
           [19.996, 24.962], [24.962, 30.026], [30.026, 39.977], [39.977, 50.065] ]

    ring_areas = [round((3.14159 * (end ** 2)) - (3.14159 * (start ** 2)), 3) for start, end in minmax]
    ring_areas = ring_areas[chstart:]
    return ring_areas


def peak_normalize_data(plotting_data, gaussian, popt):
    peak_value = gaussian(popt[1], *popt)
    return plotting_data / peak_value

def sum_list_values(lst):
    total = 0
    for value in lst:
        total += value
    return total

def half_gaussian(x, amplitude, sigma, m=1.8):
    return amplitude * np.exp(-(x - m) ** 2 / (2 * sigma ** 2)) * (x >= m)

def plot(data_dict, fiducialize_num_start=0, fiducialize_num_end=16, area_normalize=True):
    """Input a dictionary with pressure/field keys
    Output plots of Dt, sigma, fitted gaussians peak normalized, and fitted gaussians non peak normalized"""

    # Initialize variables
    plt.figure()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'lime', 'pink', 'silver', 'gray', 'maroon', 'olive', 'teal', 'navy']
    maxfittedGaussVal = []
    fittedgausses = []
    standardDevs = []
    pressures = []
    sigmasFitted = []
    normalized_data_list = []
    plotting_sets = []
    i = 0
    mean = 1.8


    mm = [0.54, 1.3, 2, 2.8, 3.9, 4.9, 5.9, 7.4, 9.9, 12.4, 14.9, 19.9, 24.85, 29.9, 39.8, 49.9]
    mm = mm[:fiducialize_num_end]
    mm = mm[fiducialize_num_start:]

    Vd = [52, 50, 48, 45, 43, 39, 35, 37, 40, 43, 45, 46, 47, 48, 49]

    for key, value in data_dict.items():
        plotting_data = None

        pressure = key.split('psi')[0].replace('neg', '-').replace('.os', '')
        pressure = round(float(pressure) * 51.715, 1)
        pressures.append(pressure)

        if area_normalize:
            plotting_data = areaNormalize(value, fiducialize_num_start)
        if not area_normalize:
            plotting_data = value

        plotting_sets.append(plotting_data)

        fittingmm, fittingdata = get_fitting_data(mean, mm, plotting_data)

        popt, pcov = fit_gaussian(half_gaussian, fittingmm, fittingdata)

        fitted_amplitude, fitted_sigma = popt
        sigmasFitted.append(fitted_sigma)

        gaussMM = np.linspace(fittingmm[0], fittingmm[-1], 1000)
        fitted_curve = half_gaussian(gaussMM, fitted_amplitude, fitted_sigma, m=mean)

        maxfittedGaussVal.append(np.max(fitted_curve))
        fittedgausses.append(fitted_curve)

        area, error = quad(lambda x: half_gaussian(x, *popt,mean), -np.inf, np.inf)
        plotting_data2 = np.array(fittingdata)
        normalized_data = plotting_data2 / area
        normalized_data_list.append(normalized_data)

        plt.scatter(mm, plotting_data, s=1, label='Fit ' + str(pressures[i]) + ' Torr')
        plt.plot(gaussMM, fitted_curve, color=colors[i], linewidth=0.5, label='Fit ' + str(pressure) + ' Torr')
        plt.title('Fitted half gaussians for mean ' + str(mean))
        plt.xlabel('Distance from center (mm)')
        plt.ylabel('Charge per mm^2')
        plt.legend(fontsize=6)
        plt.savefig('FittedGaussians', dpi=200)

        i += 1

    standardDevs = calculate_standard_deviations(mm, plotting_sets, mean)

    plot_sigma(pressures, standardDevs, mean)
    plot_sigma_fitted(pressures, sigmasFitted, mean)
    plot_diffusion_coefficient(pressures, Vd, standardDevs, 'DtData')
    plot_diffusion_coefficient(pressures, Vd, sigmasFitted, 'DtFitted')
    plot_normalized_data(normalized_data_list, pressures, fittingmm, 'chargeNormalized', 'Charge per mm^2 Normalized to 1')
    plot_normalized_data_fitted(normalized_data_list, pressures, fittingmm, 'chargeNormalizedFitted', 'Charge per mm^2 Normalized to 1', gaussMM, mean)
    plot_peak_normalized(fittedgausses, pressures, maxfittedGaussVal, gaussMM)

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

def fit_gaussian(half_gaussian, fittingmm, fittingdata):
    amplitude_guess = 1.0
    sigma_guess = 0.5
    initial_guess = [amplitude_guess, sigma_guess]
    popt, pcov = curve_fit(half_gaussian, fittingmm, fittingdata, p0=initial_guess)
    return popt, pcov

def calculate_standard_deviations(mm, plotting_sets, mean):
    standardDevs = []
    for i in range(0, len(plotting_sets)):
        ydata_norm = np.array(plotting_sets[i])
        bins = 300
        counts, bin_edges = np.histogram(mm, bins=bins, weights=ydata_norm)
        non_zero_bins = counts > 0
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        bin_centers = bin_centers[non_zero_bins]
        centerminusmean = ((bin_centers - mean) ** 2)
        countstimesdiff = counts[non_zero_bins] * centerminusmean
        nMinus1 = (sum_list_values(counts[non_zero_bins]))
        std_mm = np.sqrt(sum_list_values(countstimesdiff) / nMinus1)
        standardDevs.append(std_mm)
        print(std_mm)
    return standardDevs

def plot_sigma(pressures, standardDevs, mean):
    plt.figure()
    plt.scatter(pressures, standardDevs, label='500 & 50 V/cm')
    plt.axvline(x=745, color='grey', linestyle='--')
    plt.text(1250, 3, 'Right = 50 V/cm Left= 500 V/cm', bbox=dict(facecolor='white', edgecolor='white', boxstyle='round'))
    plt.xlabel('Pressure (Torr)')
    plt.ylabel('σ (mm)')
    plt.title('Sigma From Data for mean ' + str(mean) + ' at 500 & 50 V/cm')
    plt.savefig('SigmaData')
    plt.close()

def plot_sigma_fitted(pressures, sigmasFitted, mean):
    plt.figure()
    plt.scatter(pressures, sigmasFitted, label='500 & 50 V/cm')
    plt.axvline(x=745, color='grey', linestyle='--')
    plt.text(1250, 3, 'Right = 50 V/cm Left= 500 V/cm', bbox=dict(facecolor='white', edgecolor='white', boxstyle='round'))
    plt.xlabel('Pressure (Torr)')
    plt.ylabel('σ (mm)')
    plt.title('Sigma From Fitted Gaussian for mean ' + str(mean) + ' at 500 & 50 V/cm')
    plt.savefig('SigmaFitted')
    plt.close()

def plot_diffusion_coefficient(pressures, Vd, standardDevs, filename):
    plt.figure()
    Dtvals = []
    for i, std_mm in enumerate(standardDevs):
        Dt = (((std_mm * 0.1) ** 2) * Vd[i] * 0.1) / (2 * 11 * (1e-6))
        Dtvals.append(Dt)
        print(Vd[i])
        plt.savefig(str(i))
        plt.close()
    plt.scatter(pressures, Dtvals, label='500 & 50 V/cm')
    plt.text(1250, 17500, 'Right = 50 V/cm Left= 500 V/cm',
             bbox=dict(facecolor='white', edgecolor='white', boxstyle='round'))
    plt.axvline(x=745, color='grey', linestyle='--')
    plt.xlabel('Pressure (Torr)')
    plt.ylabel('Dt (cm^2/s)')
    plt.title('Dt Calculated from Data at 500 & 50 V/cm')
    plt.savefig(filename)
    plt.close()

def plot_normalized_data(normalized_data_list, pressures, fittingmm, filename, ylabel):
    plt.figure()
    for i in range(0, len(normalized_data_list)):
        plt.plot(fittingmm, normalized_data_list[i], linewidth=0.5, label='Fit ' + str(pressures[i]) + ' Torr')
    plt.legend(fontsize=8)
    plt.xlabel('Distance from center (mm)')
    plt.ylabel(ylabel)
    plt.title('Fitted Half Gaussians ' + ylabel + ' at 500 & 50 V/cm')
    plt.savefig(filename, dpi=1000)
    plt.close()

def plot_normalized_data_fitted(normalized_data_list, pressures, fittingmm, filename, ylabel, gaussMM,mean):
    plt.figure()
    amplitude_guess = 1.0
    sigma_guess = 0.5
    for i in range(0, len(normalized_data_list)):
        initial_guess = [amplitude_guess, sigma_guess]
        popt, pcov = curve_fit(half_gaussian, fittingmm, normalized_data_list[i], p0=initial_guess)
        fitted_amplitude, fitted_sigma = popt
        fitted_curve = half_gaussian(gaussMM, fitted_amplitude, fitted_sigma)
        plt.scatter(fittingmm, normalized_data_list[i], s=1, label='Fit ' + str(pressures[i]) + ' Torr')
        plt.plot(gaussMM, fitted_curve, linewidth=0.5, label='Fit ' + str(pressures[i]) + ' Torr')
    plt.legend(fontsize=4)
    plt.xlabel('Distance from center (mm)')
    plt.ylabel(ylabel)
    plt.title('Charge normalized Fitted half gaussians for mean ' + str(mean) + ' at 500 & 50 V/cm')
    plt.savefig(filename, dpi=1000)
    plt.close()

def plot_peak_normalized(fittedgausses, pressures,maxfittedGaussVal, gaussMM):
    plt.figure()
    maxval = np.max(maxfittedGaussVal)
    for i in range(0, len(fittedgausses)):
        plt.plot(gaussMM, normalize_data(fittedgausses[i], maxval), linewidth=0.5, label='Fit ' + str(pressures[i]) + ' Torr')
    plt.legend(fontsize=8)
    plt.xlabel('Distance from center (mm)')
    plt.ylabel('Charge per mm^2 Peak Normalized')
    plt.title('90% Ar, 10% CH4 at 500 & 50 V/cm')
    plt.savefig('peakNormalized', dpi=1000)
    plt.close()


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

def make_data_dictionary(filepaths, fiducialize_num_start=0, fiducialize_num_end=16, charge_convert=True, removeoutliers=False):
    """Looks through the list of filepaths and assumes the parent directories have the format "XXXvPerCm_posXXpXpsi"
    and parses the pressure from the name. Returns a dictionary of values to plot

    If charge_convert=true, make sure the vdd values are saved in the same directory
     as the run1_times.txt file, and in an excel sheet titled 'VddBeforeAfter.xlsx' with the Vdd values in sequential order
      from 1-16 in the second column"""

    data_dictionary = {}
    num=0

    #loop though every file, make a data dictionary for each with a key denoting pressure and feild
    for filepath in filepaths:
        num+=1

        # Initialize a list to hold 16 lists of number of resets per each channel
        reset_time_diffs = [[] for _ in range(16)]

        # Read in the data from the text file and populate the list
        with open(filepath, "r") as f:
            data = ast.literal_eval(f.read())
            for i, rtd_list in enumerate(data):
                if len(rtd_list) > 0:
                    reset_time_diffs[i] = [value for value in rtd_list if
                                          value >= 0.05]  # throw out values less than this because they aren't physical

        #do the quartile analysis to remove outliers
        if removeoutliers:
            fig, axs = plt.subplots(4, 4, figsize=(12, 12))
            axs = axs.ravel()
            for i in range(0,len(reset_time_diffs)):
                if len(reset_time_diffs[i]) != 0:
                    q1, q2, q3, upper_bound, lower_bound = calc_quartiles_bounds(reset_time_diffs[i])
                    reset_time_diffs[i] = [x for x in reset_time_diffs[i] if (x <= (upper_bound) and x >= lower_bound)]
                    #unhash to look at your outlier-removed data

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

            fig.savefig(str(num))

        if charge_convert:

            # get your Vdd values
            filepathVdds = filepath.replace('run1_times.txt', 'VddBeforeAfter.xlsx')

            # Read the Excel file into a DataFrame
            df = pd.read_excel(filepathVdds)

            # Extract values from the second column and convert them to a list
            vdd_list = df.iloc[:, 1].tolist()

            # convert the total rtds to total charge by doing cv*num_rtds which is also deltaQ*numResets=total charge
            # assume c=10pF, take vdd and multiply times 4, convert to v
            deltaQ_perReset_perChannel = [1.0e-11 * vdd * 4 * 1e-3 for vdd in vdd_list]

            # output a list of total charge seen per channel
            totalQ_PerChannel = []
            for i in range(0, 16):
                totalQ_PerChannel.append(len(reset_time_diffs[i]) * deltaQ_perReset_perChannel[i])

            # ignore everything past channel XX, this is the feducial volume cut:
            totalQ_PerChannel = totalQ_PerChannel[:fiducialize_num_end]
            totalQ_PerChannel = totalQ_PerChannel[fiducialize_num_start:]
            data = totalQ_PerChannel

        #Otherwise, just do the total number of rtds
        if not charge_convert:
            data = [len(rtd_list) for rtd_list in reset_time_diffs]

        #parse the pressure and feild from the header
        pressure = filepath.split('/')[-2].split('_')[-1].replace('psi', '').replace('p', '.')
        field = filepath.split('/')[-2].split('_')[0].replace('vPerCm', 'V/cm')

        #assign a key to your dictionary value data
        data_dictionary[pressure + "psi" + field] = data

    return data_dictionary

filepaths = get_filepaths('run1_times.txt')
data_dict= make_data_dictionary(filepaths,fiducialize_num_start=0, fiducialize_num_end=11, charge_convert=True, removeoutliers=False)
plot(data_dict,fiducialize_num_start=0, fiducialize_num_end=11, area_normalize=True)
