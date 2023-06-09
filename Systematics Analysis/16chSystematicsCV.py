"""prints CV based on number of resets seen in 5 minutes when inputting a known current source """
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

# Define the path to the text file
filelist = ["txtfiles3nA/ch1.txt","txtfiles3nA/ch2.txt","txtfiles3nA/ch3.txt","txtfiles3nA/ch4.txt","txtfiles3nA/ch5.txt","txtfiles3nA/ch6.txt","txtfiles3nA/ch7.txt","txtfiles3nA/ch8.txt","txtfiles3nA/ch9.txt","txtfiles3nA/ch10.txt","txtfiles3nA/ch11.txt","txtfiles3nA/ch12.txt","txtfiles3nA/ch13.txt","txtfiles3nA/ch14.txt","txtfiles3nA/ch15.txt","txtfiles3nA/ch16.txt"]


# t = lambda x: x * SAQ_DIV / ZYBO_FRQ
# Initialize a list to hold the 16 lists of reset time differences
reset_time_diffs = [[] for _ in range(16)]

avgrtd = []
cv=[]
v=[]

# Read in the data from the text file and populate the list
for filepath in filelist:
    chNum=int(str(filepath).split('ch')[-1].split('.')[0])

    with open(filepath, "r") as f:
        data = ast.literal_eval(f.read())
        for i, rtd_list in enumerate(data):
            if len(rtd_list) > 0:
                reset_time_diffs[i] = rtd_list

    # Loop over each channel
    for i, rtd_list in enumerate(reset_time_diffs):
        # Check if the list is empty
        if len(rtd_list) == 0 or (i+1) != chNum:
            continue
        else:
            CV=(3e-9*300)/len(rtd_list)
            cv.append(CV)


print('3Na CV value calculated using 3nA*t/numRTDS = (3e-9*300)/len(rtd_list) per channel:'+str(cv))
