#!/usr/bin/python

"""
Authors: Lisa Poole and James Pino

This pipeline is used to generate a xy scatter plot. It can be used to compare the differences between replicates
or alignment methods, etc.

"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import scipy.stats as stats

# General Information - modify appropriately for each experiment
directory = '/Users/lisapoole/Desktop/'
filename = 'HMCES_deseq_both_combined.csv'
axis_scale = 'log10' # valid options are 'log10', 'log2', or 'no' (no scaling)
entity_plotted = 'fold_change' # what you are measuring. Could be a fold_change or expression values or counts
x_axis_parameter = 'WT expression' # description of what is plotted on the X axis
y_axis_parameter = 'KO expression' # description of what is plotted on the Y axis

# Importing file(s) with data to compare
to_read = os.path.join(directory, filename) # combines the directory and filename to allow import file
data = pd.read_csv(to_read, comment='#', header=0) # imports file
print(data.dtypes) # tells what each column header is and what type of data is in each column



if axis_scale == 'log10':
    x_axis_label = "$log_10$ {} {}".format(entity_plotted, x_axis_parameter)
    y_axis_label = "$log_10$ {} {}".format(entity_plotted, y_axis_parameter)
elif axis_scale == 'log2':
    x_axis_label = "$log_2$ {} {}".format(entity_plotted, x_axis_parameter)
    y_axis_label = "$log_2$ {} {}".format(entity_plotted, y_axis_parameter)
elif axis_scale == 'no':
    x_axis_label = "{} {}".format(entity_plotted, x_axis_parameter)
    y_axis_label = "{} {}".format(entity_plotted, y_axis_parameter)
else:
    print("Error. Invalid axis scale selection. Aborting.")
    quit()

to_read = os.path.join(directory, filename) # combines the directory and filename to allow import file
data = pd.read_csv(to_read, comment='#', header=0) # imports file
print(data.dtypes) # tells what each column header is and what type of data is in each column

def plot_x_vs_y(x, y, savename, title=None):
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    x_line = np.arange(min(x),max(x))
    plt.plot(x_line,x_line,'-k', linewidth=1.5, label="ideal") # creates a line to indicate the ideal if R^2=1
    plt.plot(x, y, 'o', ms=2, color="red", label="$R^2$={0:.2f}".format(r_value**2)) # aesthetic information about plot

    if axis_scale == 'log10':
        plt.xscale('log', base=10)
        plt.yscale('log', base=10)

    elif axis_scale == 'log2':
        plt.xscale('log', base=2)
        plt.yscale('log', base=2)

    plt.xlabel(x_axis_label, fontsize=16)
    plt.ylabel(y_axis_label, fontsize=16)
    plt.title(title, fontsize=20)
    plt.annotate( "$R^2$={0:.2f}".format(r_value**2),(2,6), fontsize=16) # Include legend with R^2 correlation value
    plt.legend(loc=0, fontsize=12)
    plt.savefig("plot_{}.png".format(savename),bbox_inches='tight')
    plt.savefig("plot_{}.pdf".format(savename), bbox_inches='tight')
    plt.close()


plot_x_vs_y(data['log2FoldChange_y'], data['log2FoldChange_x'], savename='HMCES_fold_change_check',title="Comparing with and without screen control")
