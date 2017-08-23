import pandas as pd
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
# from magine.magine_analysis import Analyzer
# from magine.data.datatypes import ExperimentalData
import magine.plotting.volcano_plots as vp
from statsmodels.sandbox.stats.multicomp import multipletests

font = {'fontname':'Times New Roman'}

# data = pd.read_table('/Users/lisapoole/Desktop/E65/E65_sigchange.txt')
data = pd.read_table('/Users/lisapoole/Desktop/HMCES_SE_screen_sig_changed.txt')
data['data_type'] = 'rna_seq'

fold_change = 'treated_control_fold_change'
# data[fold_change] = data['logFC']
data[fold_change] = data['log2FoldChange']
flag = 'significant_flag'

# data = data[data['padj'].notnull()]
# # data['p_value_group_1_and_group_2'] = data['FDR']
# data['p_value_group_1_and_group_2'] = data['padj']
# pvals = data['padj']
# print(np.shape(pvals))
# adj_pvalues = multipletests(pvals=pvals, method='fdr_bh')[1]
# print(np.shape(adj_pvalues))
# print(adj_pvalues)
# data['p_value_group_1_and_group_2'] = adj_pvalues
# import matplotlib.pyplot as plt
# plt.plot(data['padj'], data['p_value_group_1_and_group_2'], '.')
# plt.show()

data['p_value_group_1_and_group_2'] = data['padj']
# data['p_value_group_1_and_group_2'] = data['PValue']
data[fold_change] = np.where(data[fold_change] > 2,
                             np.exp2(data[fold_change]), -1*np.exp2(-1*data[fold_change]))
data[flag] = False

location = (np.abs(data[fold_change]) > 1.5) & (data['padj'] < 0.1)
data.loc[location, flag] = True
data['species_type'] = 'protein'
data['gene'] = data['gene_id']
data['protein'] = data['gene']
data['time'] = 0

data.to_csv('HMCES_two_controls_data.csv', index=False)
# print(data.dtypes)
# exp_data = ExperimentalData('magine_RPE_data.csv', data_directory='.')
vp.volcano_plot(data, 'HMCES_controls', image_type='pdf')#, x_range=[-3,3], y_range=[-1,12])

quit()

magine = Analyzer(exp_data,network='1',metric='enrichment', output_directory='RPE_HTML')
magine.run_go(data_type='rnaseq', run_type='changed', metric='enrichment')
quit