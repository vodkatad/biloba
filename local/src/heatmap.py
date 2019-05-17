import sys, os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.colors as colors

if len(sys.argv) < 3:
    print("Missing Arguments")
    exit(1)
    
if os.path.isfile(sys.argv[1]) == False:
    print("FileNotFoundException")
    exit(1)

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)
    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

cnvs = pd.read_csv(sys.argv[1], sep="\t", index_col=0, usecols = lambda column : column not in ['START', 'END'])
cnvs = cnvs.transpose()

norm=MidpointNormalize(vmin=0., midpoint=2., vmax=12.)
cbar_kws={"ticks":np.arange(0,13,1)}

#hm = sns.clustermap(cnvs, cmap='RdBu_r', robust=True, col_cluster=False, yticklabels=False, annot_kws={'size':12}, vmax=16, center=2)
sns.clustermap(cnvs, method='complete', metric='euclidean', 
        col_cluster=False, yticklabels=False, cmap='RdBu_r', 
        vmin=0, vmax=12,
        norm=norm,
        cbar_kws=cbar_kws)

plt.gcf().set_size_inches(18.5, 10.5)
plt.savefig(sys.argv[2])
