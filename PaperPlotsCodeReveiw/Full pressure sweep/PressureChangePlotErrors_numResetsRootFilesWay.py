import matplotlib.pyplot as plt
import ast
import numpy as np
import os
import pandas as pd
from array import array
import ROOT
from ROOT import TGraph
from ROOT import TCanvas, TPaveText, TPad,kRed, TF2, TGraph, TGraphErrors, TF1, gRandom, TVirtualFitter
from ROOT import gROOT, gStyle
#from scipy import signal

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


def weighted_average_time_bins(time_list, num_bins=50):
    """Calculates a weighted average based on the number of times in each bin, using a default 50 bins"""
    if len(time_list) == 0:
        weighted_mean = 0
    else:
        # Calculate the bin edges
        bin_edges = np.linspace(min(time_list), max(time_list), num_bins + 1)

        # Bin the values and get the bin counts
        binned_values, _ = np.histogram(time_list, bins=bin_edges)

        # Get the bin indices for each value in time_list
        bin_indices = np.digitize(time_list, bin_edges) - 1

        # Multiply each value by its weight
        weighted_values = [time_list[i] * binned_values[bin_indices[i]] if bin_indices[i] < num_bins else 0 for i in
                           range(len(time_list))]

        # Calculate the weighted sum of 'weighted_values'
        weighted_sum = sum(weighted_values)

        # Calculate the sum of weights, which involves summing values from 'binned_values' based on 'bin_indices'
        # If 'bin_indices' is out of bounds (greater than or equal to 'num_bins'), use 0 for that index
        sum_weights = sum(
            [binned_values[bin_indices[i]] if bin_indices[i] < num_bins else 0 for i in range(len(time_list))])

        # Calculate the weighted mean by dividing the weighted sum by the sum of weights
        weighted_mean = weighted_sum / sum_weights

    return weighted_mean


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
            areaErrFinal=[charge * areaer for charge, areaer in zip(plotting_data,areaErr)]

            # Add the area error above and below each data point
            yerrAbove = [yerrAb + areaEr for yerrAb,areaEr in zip(yerrAbove,areaErrFinal)]
            yerrBelow = [yerrbel + areaEr for yerrbel,areaEr in zip(yerrBelow,areaErrFinal)]

            if plot_charge:
                ax.set_ylabel('Total Charge (C) per mm\u00b2')
            if not plot_charge:
                ax.set_ylabel('Number of times per mm^2')

        if not area_normalize:

            plotting_data = value

            if plot_charge:
                ax.set_ylabel('Total Charge (C)')
            if not plot_charge:
                ax.set_ylabel('Number of times')

        # apply your filters if you only want to plot above 1 atm:
        if round(float(pressure) * 51.715, 1) > 720:
            ax.plot(mm, plotting_data, marker='.', label=str(round(float(pressure) * 51.715, 1)) + ' Torr', color=colors[n])
            ax.errorbar(mm, plotting_data, yerr=(yerrBelow,yerrAbove), fmt='none', capsize=6, markersize='6', color=colors[n])
        n += 1

    ax.set_xlabel('Radius (mm)')
    ax.set_title(title)
    ax.legend()
    fig.savefig(save_file_name)

def make_rootPlot(xTimeData,yResetsData, pressure):
    c1 = TCanvas('c1', 'Example of plot', 200, 10, 700, 900)
    mg = ROOT.TMultiGraph('mg', 'TMultiGraph')

    markernums=[34,29,21,22,23,24,25,26,27,28,29,30,31,32,33,34]
    colors=[1,2,3,4,5,6,7,8,9,38,36,46,32,32,40,38,28,12,48]
    for number in range(0,len(xTimeData)):
        if len(xTimeData[number]) == 0:
            continue
        else:
            graph = TGraph(len(xTimeData[number]),array('d',xTimeData[number]),array('d',yResetsData[number]))
            graph.SetName("Channel "+str(number+1));
            graph.SetTitle("Channel "+str(number+1));
            graph.SetMarkerStyle(markernums[number])
            graph.SetMarkerSize(1)
            graph.SetMarkerColor(colors[number])
            graph.Fit("pol1", "FQ")
            graph.SetLineColor(colors[number])
            graph.GetFunction("pol1").SetLineColor(colors[number])

            graph.Draw("AP")
            mg.Add(graph)
    mg.SetTitle(str(round(float(pressure), 1)) + " Torr")
    mg.Draw("AP")
    mg.GetXaxis().SetTitle("Time (s)")
    mg.GetYaxis().SetTitle("Total Number of Resets")
    legend = c1.BuildLegend(0.1, 0.7, 0.4, 0.9)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.02)
    c1.SaveAs(str(round(float(pressure), 1)).replace('.','p') + ".png")

def make_rootMeanAndErrs(TimeData, ResetsData, pressure):
    '''for i in range(0,len(TimeData)):
        listOfResetTimeStamps = TimeData[i]
        rtds=[]
        timestamp=0
        for n in listOfResetTimeStamps:
            timeDiff = n-timestamp
            rtds.append(timeDiff)
            timestamp=n'''

    #sos = signal.butter(10, 1 / 50, fs=1 / 50 * 3, output='sos')
    #filtered = signal.sosfilt(sos, timestamp_data)

    c1 = TCanvas('c1', 'Example of plot', 200, 10, 700, 900)
    mg = ROOT.TMultiGraph('mg', 'TMultiGraph')
    #This splits plots into subplots if you want to plot everything on same plot do c1.cd(1,1)
    c1.Divide(4, 4)


    markernums = [34, 29, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34]
    colors = [1, 2, 3, 4, 5, 6, 7, 8, 9, 38, 36, 46, 32, 32, 40, 38, 28, 12, 48]

    y_position = 400  # Initial y-position

    # go through each channel list. the passed x and y data have 16 lists, one for each channel, for one pressure data
    for number in range(0, len(TimeData)):

        # only plot the channels with data
        if len(TimeData[number]) == 0:
            continue
        else:
            graph = TGraph(len(TimeData[number]), array('d', ResetsData[number]), array('d', TimeData[number]))
            graph.SetName("Channel " + str(number + 1))
            graph.SetTitle("Channel " + str(number + 1))
            graph.SetMarkerStyle(markernums[number])
            graph.SetMarkerSize(1)
            graph.SetMarkerColor(colors[number])

            stopFitIndex=[]
            startNextFitIndex=[]
            timestamp = 0
            for n in range(0, len(TimeData[number])):
                if n == 0:
                    timestamp = TimeData[number][n]

                timeDiff = TimeData[number][n] - timestamp
                timestamp = TimeData[number][n]
                totDiff = TimeData[number][-1] - TimeData[number][0]
                metric = totDiff/5
                if timeDiff>metric:
                    stopFitIndex.append(n)
                    startNextFitIndex.append(n+1)
            if len(startNextFitIndex) == 0:
                graph.Fit("pol1", "FQ")
            else:
                for index in range(0, len(stopFitIndex)):
                    if index == 0:
                        graph.Fit("pol1", "FQ", "", 0, ResetsData[number][stopFitIndex[index]])
                    if index == stopFitIndex[-1]:
                        graph.Fit("pol1", "FQ", "", ResetsData[number][startNextFitIndex[index]],
                                  ResetsData[number][-1])
                    else:
                        graph.Fit("pol1", "FQ", "", ResetsData[number][startNextFitIndex[index]], ResetsData[number][stopFitIndex[index+1]])

            #if (number+1) == 9 or (number+1) == 6:
            #    graph.Fit("pol1", "FQ", "", 20, ResetsData[number][-1])
            #else:
            '''
            d = np.diff(TimeData[number])
            b=d>50
            inds=np.where(b)
            print(inds)

            f1 = ROOT.TF1("f", "pol1", ResetsData[number][0], inds[0])
            graph.Fit(f1, "FR")
            f1.Draw("SAME")
            graph.Fit("pol1","","", inds[0]+2, ResetsData[number][-1])
            '''
            #graph.Fit("pol1", "FQ")


            graph.SetLineColor(colors[number])
            graph.GetFunction("pol1").SetLineColor(colors[number])

            fit_result = graph.GetFunction("pol1")
            # this is 1/slope where slope is the value of avg rtd
            slope = fit_result.GetParameter(1)

            # Create a TLatex object and set its position relative to the graph
            latex = ROOT.TLatex(80, y_position, f"Channel {number + 1}: Slope: {slope:.4f}")
            latex.SetTextSize(0.02)
            latex.SetTextColor(colors[number])
            graph.GetListOfFunctions().Add(latex)  # Add the latex label to the graph

            y_position += 50 # Increment the y-position for the next label

            #make a araay of points for error bars thats the same length as x vals
            ngr=len(ResetsData[number])
            errorPoints = TGraphErrors(ngr)
            for i in range(ngr):
                xvals = graph.GetX()[i]
                errorPoints.SetPoint(i, xvals , 0);

            TVirtualFitter.GetFitter().GetConfidenceIntervals(errorPoints)
            errorPoints.SetLineColor(kRed);
            errorPoints.Draw("AP")

            graph.SetMarkerStyle(5)
            graph.SetMarkerSize(0.7)
            graph.Draw("psame")
            mg.Add(graph)

    mg.SetTitle(str(round(float(pressure), 1)) + " Torr")
    mg.Draw("AP")
    mg.GetYaxis().SetTitle("Time (s)")
    mg.GetXaxis().SetTitle("Total Number of Resets")
    legend = c1.BuildLegend(0.1, 0.7, 0.4, 0.9)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.02)

    c1.SaveAs("Filtered" + str(round(float(pressure), 1)).replace('.', 'p') + ".png")



def get_filepaths(textfile, current_directory=True):
    """input a text file containing a 2d list with 16 lists of times, parsed from
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

                if "pycache" in str(dir_path):
                       continue
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

        # Initialize a list to hold the 16 lists of total number of resets at each time (will be 12345...end)
        totNumResets = [[] for _ in range(16)]

        # Initialize a list to hold the 16 lists of time of each reset
        timeOfResets = [[] for _ in range(16)]

        # Initialize a list to hold the 16 lists of reset time differences for backgrounds
        background_reset_time_diffs = [[] for _ in range(16)]
        background_data = None

        # Get your vdd file names
        filepathVdds = filepath.replace('run1_times.txt', 'VddBeforeAfter.xlsx')

        # Read the Excel file into a DataFrames
        df = pd.read_excel(filepathVdds)

        # Extract values from the second column and convert them to a list
        vdd_list = df.iloc[:, 1].tolist()

        # Read in the data from the text file and populate the list
        with open(filepath, "r") as f:
            # Read the contents of the file and use 'ast.literal_eval()' to safely parse it as a Python literal
            data = ast.literal_eval(f.read())

            # Iterate over the elements of 'data' using enumeration
            num2=0
            for i, time_list in enumerate(data):
                num2+=1
                vddVal = vdd_list[i]
                # Check if there are times, if not leave the list empty as it already is
                #time_list = [value for value in time_list if value >= 0.05 * vddVal * 1e-3 * 4]
                if len(time_list) > 0:
                    # Make a list of the total number of resets (just a list from 1 to len(resets)
                    #totNumResets[i] = list(range(1, len(time_list) + 1))


                    indexTimestamp=0
                    indexTime = 0
                    numResets=0
                    index=0
                    for timestamp in time_list:
                        time = timestamp * 200 / 30e6
                        timeDiff=time-indexTime
                        if timeDiff <50e-3: #throw out nonphysical errors due to sparking
                            continue
                        if indexTimestamp == timestamp:
                            continue
                        else:
                            index+=1
                            numResets+=1
                            indexTimestamp = timestamp
                            indexTime = time
                            if index<5:
                                continue
                            else:
                                totNumResets[i].append(numResets)
                                timeOfResets[i].append(time)  # Append this time
                else:
                    continue
        pressure = filepath.split('/')[-2].split('_pos')[-1].replace('psi', '').replace('p', '.')
        make_rootPlot(timeOfResets, totNumResets, str(float(pressure) * 51.715))
        make_rootMeanAndErrs(timeOfResets, totNumResets, str(float(pressure) * 51.715))

        # Define the number of rows and columns for subplots
        num_rows, num_cols = 4, 4

        # Create the subplots with tight layout
        fig, axs = plt.subplots(num_rows, num_cols, figsize=(12, 12), constrained_layout=True)

        # Flatten the axs array for easier iteration
        axs = axs.ravel()

        # Loop through your data
        for i in range(0, len(timeOfResets)):
            # Check if there is more than one reset data point to avoid division by 0 error
            if len(timeOfResets[i]) > 1:
                axs[i].scatter(timeOfResets[i], totNumResets[i])
                axs[i].set_xlabel("Time (s)")
                axs[i].set_ylabel("Total number of resets")
                #axs[i].legend(fontsize=6)
                axs[i].set_title("Ch{}".format(i + 1))
                #if (i + 1)==7 and num==4:
                #    print(data[i])

        fig.suptitle('Pressure: '+str(float(pressure)* 51.715) + ' Torr', fontsize=16)

        # Automatically adjust subplot spacing to prevent overlap
        plt.tight_layout()
        #plt.show()

        # Save the figure
        fig.savefig(str(num))


filepaths = get_filepaths('run1_times.txt')
make_data_dictionary(filepaths, fiducialize_num_start=0,fiducialize_num_end=11, background_subtracted=False,charge_convert=True)
