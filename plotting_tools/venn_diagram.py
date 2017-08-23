import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2 as _venn2

directory = '/Users/lisapoole/Desktop'
number_of_comparisons = 2 # valid options are 1 or 2
filename1 = 'VENN_END.txt'
filename2 = 'VENN_END1.txt'
column1 = 'snp' # column name of columns used for analysis
column2 = 'snp' # column name of columns used for analysis
column_number1 = 0 # Column numbers from file used for analysis, start numbering with 0
column_number2 = 0 # Column numbers from file used for analysis, start numbering with 0

if number_of_comparisons ==2:
    # Import data
    to_read1 = os.path.join(directory, filename1) # combines the directory and filename to allow import file
    to_read2 = os.path.join(directory, filename2) # combines the directory and filename to allow import file
    data1 = pd.read_table(to_read1, comment='#', usecols=[column_number1],names=[column])
    data2 = pd.read_table(to_read2, comment='#', usecols=[column_number2],names=[column])

def create_venn2(list1, list2, label1, label2, save_name,
                 out_dir=None, image_format='svg'):
    """

    Parameters
    ----------
    list1 : list_like
    list2 : list_like
    label1 : str
    label2 : str
    save_name : str
    out_dir : str, optional
    image_format : str, optional
        default png

    Returns
    -------

    """
    set1 = set(list1)
    set2 = set(list2)
    v = _venn2([set1, set2], ('', ''), #('%s(%s)' % (label1, str(len(set1))),
                              #'%s(%s)' % (label2, str(len(set2))),)
               set_colors=('blue','yellow')
               )
    plt.tight_layout()
    save_name = "{}.{}".format(save_name, image_format)
    if out_dir is not None:
        save_name = os.path.join(out_dir, save_name)
    plt.savefig(save_name, bbox_inches='tight')
    plt.close()

create_venn2(list(sample_1), list(sample_2),'', '', 'sm1_with_without_hu')
set_1 = set(sample_1)
set_2 = set(sample_2)

set_2 = set_1.intersection(set_2)
file_output = ''
for i in sorted(set_2):
    print(i)
    file_output+=i+'\n'

with open('sm1_with_without_treat.txt', 'w') as f:
    f.write(file_output)