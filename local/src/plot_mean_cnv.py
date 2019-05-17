import sys, os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

if len(sys.argv) < 3:
    print("Missing Arguments")
    exit(1)
if os.path.isfile(sys.argv[1]) == False:
    print("FileNotFoundException")
    exit(1)

cnvs = pd.read_csv(sys.argv[1], sep="\t", index_col=0, usecols = lambda column : column not in ['START', 'END'])
                       
average = cnvs.mean(axis=1)
cnvs['mean_cnv'] = average.values
cnvs.plot(y='mean_cnv', grid=True, title = 'Mean cnv profile', legend=False)

plt.gcf().set_size_inches(18.5, 10.5)
plt.savefig(sys.argv[2]dpi=100)

