"""Saves a table for various pressures and feilds on the drift with the flashlamp off"""
from tabulate import tabulate
import numpy as np

def calc_quartiles_bounds(data):
    q1 = np.percentile(data, 25)
    q2 = np.percentile(data, 50)
    q3 = np.percentile(data, 75)

    iqr = q3 - q1
    upper_bound = q3 + 1.5 * iqr
    lower_bound = q1 - 1.5 * iqr

    return q1, q2, q3, upper_bound, lower_bound
values_list25 = []
currentFiles25=['Picoameter/background/25v_per_cm/26.7psi/current_time0.txt','Picoameter/background/25v_per_cm/27.7psi/current_time0.txt','Picoameter/background/25v_per_cm/28.7psi/current_time0.txt','Picoameter/background/25v_per_cm/29.7psi/current_time0.txt']
values_list25std = []
for filepath in currentFiles25:
    values_list=[]
    with open(filepath, "r") as f:
        for line in f:
            # Split the line by space and take the first column
            columns = line.split()
            if columns:
                value = columns[0]
                values_list.append(float(value))
    q1, q2, q3, upper_bound, lower_bound = calc_quartiles_bounds(values_list)
    values_list = [x for x in values_list if (x <= (upper_bound) and x >= lower_bound)]
    values_list25.append(np.mean(values_list))
    values_list25std.append(np.std(values_list))

#-----------------------------------------------------------------------------------------
values_list50 = []
currentFiles50=['Picoameter/background/50v_per_cm/26.7psi/current_time0.txt','Picoameter/background/50v_per_cm/27.7psi/current_time0.txt','Picoameter/background/50v_per_cm/28.7psi/current_time0.txt','Picoameter/background/50v_per_cm/29.7psi/current_time0.txt']
values_list50std = []
for filepath in currentFiles50:
    values_list=[]
    with open(filepath, "r") as f:
        for line in f:
            # Split the line by space and take the first column
            columns = line.split()
            if columns:
                value = columns[0]
                values_list.append(float(value))
    q1, q2, q3, upper_bound, lower_bound = calc_quartiles_bounds(values_list)
    values_list = [x for x in values_list if (x <= (upper_bound) and x >= lower_bound)]
    values_list50.append(np.mean(values_list))
    values_list50std.append(np.std(values_list))
#-----------------------------------------------------------------------------------------------------
values_list60 = []
currentFiles60=['Picoameter/background/60v_per_cm/26.7psi/current_time0.txt','Picoameter/background/60v_per_cm/27.7psi/current_time0.txt','Picoameter/background/60v_per_cm/28.7psi/current_time0.txt','Picoameter/background/60v_per_cm/29.7psi/current_time0.txt']
values_list60std = []
for filepath in currentFiles60:
    values_list=[]
    with open(filepath, "r") as f:
        for line in f:
            # Split the line by space and take the first column
            columns = line.split()
            if columns:
                value = columns[0]
                values_list.append(float(value))
    q1, q2, q3, upper_bound, lower_bound = calc_quartiles_bounds(values_list)
    values_list = [x for x in values_list if (x <= (upper_bound) and x >= lower_bound)]
    values_list60.append(np.mean(values_list))
    values_list60std.append(np.std(values_list))

# Data for the table
data = [
    ["25 V/cm", str(round(values_list25[0],16)) + str(' \u00B1 ') + str(round(values_list25std[0],16)) +' A', str(round(values_list25[1],16))+str(' \u00B1 ') + str(round(values_list25std[1],16))+ ' A', str(round(values_list25[2],16))+str(' \u00B1 ') + str(round(values_list25std[2],16))+ ' A', str(round(values_list25[3],16))+str(' \u00B1 ') + str(round(values_list25std[3],16))+ ' A'],
    ["50 V/cm",str(round(values_list50[0],16))+str(' \u00B1 ') + str(round(values_list50std[0],16))+ ' A', str(round(values_list50[1],16))+str(' \u00B1 ') + str(round(values_list50std[1],16))+ ' A', str(round(values_list50[2],16))+str(' \u00B1 ') + str(round(values_list50std[2],16))+ ' A', str(round(values_list50[3],16))+str(' \u00B1 ') + str(round(values_list50std[3],16))+ ' A'],
    ["60 V/cm",str(round(values_list60[0],16))+str(' \u00B1 ') + str(round(values_list60std[0],16))+ ' A', str(round(values_list60[1],16))+str(' \u00B1 ') + str(round(values_list60std[1],16))+ ' A', str(round(values_list60[2],16))+str(' \u00B1 ') + str(round(values_list60std[2],16))+ ' A', str(round(values_list60[3],16))+str(' \u00B1 ') + str(round(values_list60std[3],16))+ ' A']
]

# Table headers
headers = ["", "26.7 PSI", "27.7 PSI", "28.7 PSI","29.7 PSI"]

# Title of the table
title = "Background Electric Field vs. Pressure (68% CI)"

# Generate the table with title
table = tabulate(data, headers, tablefmt="rounded_grid")

# Add the title at the top of the table
table_with_title = f"{title}\n{table}"

# Save the table to a file
with open("table.txt", "w") as file:
    file.write(table_with_title)

