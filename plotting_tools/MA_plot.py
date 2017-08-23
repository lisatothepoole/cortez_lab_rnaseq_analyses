#!/usr/bin/python

""" Create an MA plot of RNAseq data
Authors: Lisa Poole and James Pino

This pipeline is optimized for output from DeSeq, which already has the baseMean calculated.

"""

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd


def ma_plot(data, savename, out_dir=None, save_format='png', x_range=None, y_range=None, significance_parameter='padj',
            p_value_cutoff=0.05):
    """
    Creates an MA plot of data type provided

    Parameters:
    ___________

    data: pandas.Dataframe
        data to create MA plot from
    savename: str
        name to save figure to
    out_dir: str, directory, optional
        location to save figure
    save_format: str
        format of plot to save {'png,' 'pdf,' 'svg,' etc}
    x_range: array_like, optional
        upper and lower bounds of plot in x direction
    y_range: array_like, optional
        upper and lower bounds of plot in y direction
    significance_parameter: str, optional
        significance parameter to plot {adjusted p value 'padj' or p value 'pval'}
    p_value_cutoff: float, optional
        Criteria for significance

    """

    data = data[data['foldChange'].notnull()]
    data['log_mean_values'] = np.log10(data['baseMean'])

    crit_1 = data[significance_parameter] < p_value_cutoff
    crit = (crit_1)

    non_crit = data[~crit]  # makes a list of genes that aren't significantly changed
    crit = data[crit]  # makes a list of genes that are significantly changed

    fig = plt.plot(non_crit['log_mean_values'], non_crit['log2FoldChange'], 'o')
    fig = plt.plot(crit['log_mean_values'], crit['log2FoldChange'],
                   'or')  # plot significantly changed values as a red dot
    ax = plt.subplots()
    if y_range is not None:
        ax.set_ylim(y_range[0], y_range[1])
    if x_range is not None:
        ax.set_xlim(x_range[0], x_range[1])
    plt.ylabel('$log_2$(Fold Change)')
    plt.xlabel('$log_10$(mean_counts)')

    save_plot(fig, save_name=savename, out_dir=out_dir, image_type=save_format)


def save_plot(fig, save_name, out_dir=None, image_type='png'):
    """

    Parameters
    ----------
    fig : matplotlib.pyplot.figure
        Figure to be saved
    save_name : str
        output file name
    out_dir : str, optional
        output path
    image_type : str, optional
        output type of file, {"png", "pdf", etc..}

    Returns
    -------

    """
    fig.tight_layout()
    save_name = '{}.{}'.format(save_name, image_type)
    if out_dir is not None:
        save_name = os.path.join(out_dir, save_name)
    fig.savefig(save_name, bbox_inches='tight')
    plt.close()
