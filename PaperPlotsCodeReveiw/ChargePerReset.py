import matplotlib.pyplot as plt
import ast
import numpy as np
import pandas as pd

def calc_quartiles_bounds(data):
    q1 = np.percentile(data, 25)
    q2 = np.percentile(data, 50)
    q3 = np.percentile(data, 75)

    iqr = q3 - q1
    upper_bound = q3 + 1.5 * iqr
    lower_bound = q1 - 1.5 * iqr

    return q1, q2, q3, upper_bound, lower_bound


# Define the path to the text file
filelist1 = ["2nA/2nA/txtfiles2nA/ch1.txt","2nA/2nA/txtfiles2nA/ch2.txt","2nA/2nA/txtfiles2nA/ch3.txt","2nA/2nA/txtfiles2nA/ch4.txt","2nA/2nA/txtfiles2nA/ch5.txt","2nA/2nA/txtfiles2nA/ch6.txt","2nA/2nA/txtfiles2nA/ch7.txt","2nA/2nA/txtfiles2nA/ch8.txt","2nA/2nA/txtfiles2nA/ch9.txt","2nA/2nA/txtfiles2nA/ch10.txt","2nA/2nA/txtfiles2nA/ch11.txt","2nA/2nA/txtfiles2nA/ch12.txt","2nA/2nA/txtfiles2nA/ch13.txt","2nA/2nA/txtfiles2nA/ch14.txt","2nA/2nA/txtfiles2nA/ch15.txt","2nA/2nA/txtfiles2nA/ch16.txt"]
currentFiles= ['drive-download-20230604T191134Z-001/2nA/currentSourceUncertainty/current_time0.txt']

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

# Initialize a list to hold the 16 lists of reset time differences
reset_time_diffs = [[] for _ in range(16)]

numRtds2nA = []
absoluteUncertainty2nA=0
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
            numRtds2nA.append(len(rtd_list))
            absoluteUncertainty2nA=np.sqrt(len(rtd_list))


fig, axs = plt.subplots(4, 4, figsize=(12, 12))
axs = axs.ravel()
excelname = '2nA/VddsBeforeAfterCalibrationRun.xlsx'

# Read the Excel file into a DataFrame
df = pd.read_excel(excelname)

# Extract values from the second column and convert them to a list
vddch = df.iloc[:, 1].tolist()

QperReset=[]
for i in range(0, 16):
    xnumRtds = [numRtds2nA[i]]
    xer = [[absoluteUncertainty2nA], [absoluteUncertainty2nA]]
    y = 1e-9*300
    QperReset.append(y/xnumRtds[0])



plt.figure()
ch=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
ch=ch[:11]
QperReset=QperReset[:11]
plt.title('Average Charge Per Reset for 2nA input current over 5 minutes')
plt.ylabel('Q/Reset (C/Reset)')
plt.xlabel('Channel Number')
plt.plot(ch,QperReset, label = 'Charge/Reset')
plt.legend(loc='upper left')

plt.twinx()

# Plot the second dataset on the second y-axis
vddch = vddch[:11]
plt.plot(ch, vddch, color = 'green', label = 'Vdd Value')
plt.ylabel('Vdd (mV)')
plt.legend(loc='upper right')

#plt.errorbar(xnumRtds, y, yerr=yer, xerr=xer, fmt='o', markersize=2, linewidth=0.01, capsize=8, color='black', ecolor='red')
plt.savefig('ChargePerReset')
