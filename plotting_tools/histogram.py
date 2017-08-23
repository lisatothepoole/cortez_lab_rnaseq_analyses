#!/usr/bin/python

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

# General Information - modify appropriately for each experiment
directory = '/Users/lisapoole/Desktop/'
filename = 'Cortez-U2OS-siNT_blacklist_filtered_peaks.csv'
histogram_cutoff = 'none' # can select a number that everything above it is binned together. If no cutoff is desired, insert 'none'
histogram_cutoff_add_1 = 101 # only applicable if number is selected in histogram_cutoff
bin_sizes = 'auto' # if you prefer auto, insert 'auto'
half_bin_size = 5 # only valid if bin_sizes is a number
column_header = 'pileup' # name of the column pulled from file for histogram
log_axis = 'no' # valid options are 'yes' for log transformation of y axis or 'no' for no log transformation
axis_scale = 'none' # valid options are 10 for log10, 2 for log2, or 'none'
save_format = 'pdf' # file extension for saving files. common options are 'pdf,' '.svg,' or '.png'

# Importing file(s) with data to compare
to_read = os.path.join(directory, filename) # combines the directory and filename to allow import file
data = pd.read_csv(to_read, comment='#', header=0) # imports file
print(data.dtypes) # tells what each column headers and types of data each column is
data[column_header] = data[column_header].apply(pd.to_numeric) # converts the column to numbers in case they are not


def peak_histogram_plot(savename, title):

    if bin_sizes == 'auto':
        bins = 'auto'
    else:
        bins = range(0, histogram_cutoff_add_1, bin_sizes)


    if histogram_cutoff == 'none':
        data_for_hist = data.copy()
        fig, ax = plt.subplots(figsize=(9, 5))
        plt.hist(data_for_hist[column_header],
                 color='red', bins=bins)

        bin_labels = np.array(bins, dtype='|S4')
        print(bin_labels)
        # bins = np.array(bins) + 5
        # plt.xticks(bins)
        # bins = np.array(bins) + 5
        print(bins)
        ax.set_xticklabels(bins, fontsize=10)

    else:
        data_for_hist = data.copy()
        cutoff = histogram_cutoff
        data_for_hist.loc[data_for_hist[column_header] > cutoff, column_header] = cutoff + 1

        bins.append(cutoff + bin_sizes)
        print(bins)
        fig, ax = plt.subplots(figsize=(9, 5))
        plt.hist(data_for_hist[column_header], bins=bins,
                           color='red')
        bin_labels = np.array(bins, dtype='|S4')
        bin_labels = []
        for i in bins:
            if i == histogram_cutoff:
                bin_labels.append('>100')
            elif i >histogram_cutoff:
                continue
            else:
                bin_labels.append('{}-{}'.format(i,i+bin_sizes))
        print(bin_labels)

        bins = np.array(bins)+half_bin_size
        print(bins)

        ax.set_xticklabels(bin_labels, fontsize=10)


    if log_axis == 'yes':
        plt.yscale('log', base=axis_scale)
        plt.ylabel('$log_{}$ Frequency'.format(axis_scale))
    elif log_axis == 'no':
        plt.ylabel('Frequency'.format(axis_scale))


    if bin_sizes == 'auto':
        pass
    else:
        ax.set_xticks(bins)

    # ax.set_xticks(bins)
    plt.xlim(0, 115)
    plt.title('{}'.format(title), fontsize=20)
    plt.xlabel('{}'.format(column_header))
    plt.savefig("{}/{}.{}".format(directory, savename, save_format), bbox_inches='tight')
    # plt.show()


peak_histogram_plot('practice_histogram', 'practice_histogram')

