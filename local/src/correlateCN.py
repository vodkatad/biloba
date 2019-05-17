import sys
import argparse
import numpy as np
from scipy.stats import spearmanr
import h5py as h5
import pandas as pd
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description="Spearman correlation on a given list of cells")

parser.add_argument("lower", metavar='L', action='store',
                        help='lower threshold of the set of ploidies to discard',
                                                nargs=1, type=float)

parser.add_argument("higher", metavar='H', action='store',
                        help='higher threshold of the set of ploidies to discard',
                                                nargs=1, type=float)
parser.add_argument("filter", choices=['tumoral','normal'], metavar='tumoral/normal', action='store',
                        help='The cell subset to be taken in account',
                        nargs=1, type=str)

parser.add_argument("hdf5", metavar='cnv_data.h5', action='store',
                help='hdf5 file storing 10x intermediate results',
                        nargs=1, type=str)

parser.add_argument("cells_summary", metavar='per_cell_summary_metrics.csv', action='store',
                        help='csv results storing 10x per cell summary results',
                                                nargs=1, type=str)

parser.add_argument("segcopy", metavar='SegCopy', action='store',
                help='csv file storing Ginkgo binned results',
                        nargs=1, type=str)

parser.add_argument("results", metavar='results.txt', action='store',
                        help='txt file storing Ginkgo per cell summary results',
                                                nargs=1, type=str)

parser.add_argument("sample", metavar='Sample1', action='store',
                                help='Sample name',
                               nargs=1, type=str)

parser.add_argument("outdir", metavar='/path/to/out/dir', action='store',
                help='path to the desired output directory', 
                        nargs=1, type=str)

args=parser.parse_args()

h5_f = args.hdf5[0]
ginkgo_f = args.segcopy[0]
stats_f = args.cells_summary[0]
results_f = args.results[0]
sample = args.sample[0]
outdir = args.outdir[0]

low_thresh = args.lower[0]
high_thresh = args.higher[0]
filter = args.filter[0]


# select valid cells
results_df = pd.read_csv(results_f, sep="\t")
results_df = results_df.dropna()

stats_df = pd.read_csv(stats_f, usecols=['barcode', 'cell_id', 'mean_ploidy', 'is_noisy'])
stats_df['barcode'] = stats_df['barcode'].str[:-2]

cells = []
if filter == "tumoral":
    cells_g = results_df['Sample'].loc[(results_df["Copy_Number"] < low_thresh) | (results_df["Copy_Number"] > high_thresh)].values
    cells_10ex = stats_df['barcode'].loc[(stats_df["mean_ploidy"] < low_thresh) | (stats_df["mean_ploidy"] > high_thresh)].values
    cells = cells_g
    #union between the subsets of tumor cells according to both tools
    for cell in cells_10ex:
        if cell in cells:
            np.append(cells, np.array([cell]))
else:
    cells_g = results_df['Sample'].loc[(results_df["Copy_Number"] >= low_thresh) & (results_df["Copy_Number"] <= high_thresh)].values
    cells_10ex = stats_df['barcode'].loc[(stats_df["mean_ploidy"] >= low_thresh) & (stats_df["mean_ploidy"] <= high_thresh)].values
    #intersection between the  subsets of normal cells according to both tools
    for cell in cells_g:
        if cell in cells_10ex:
            cells.append(cell)

num_cells = len(cells)

with open(outdir+"/cell_stats.txt", "w+") as f:
    f.write("Ginkgo_cells_N\t10x_cells_N\tIntersection_size\n")
    f.write(str(len(cells_g)) + "\t" + str(len(cells_10ex)) + "\t" + str(num_cells) + "\n")


# Load Ginkgo data
df_ginkgo = pd.read_csv(ginkgo_f, sep="\t")
chr_list = df_ginkgo['CHR'].unique()
all_cells = df_ginkgo.columns[4:]


# Load 10x data
h5_data = h5.File(h5_f)
tot_cells = h5_data["constants/num_cells"].value
df_10ex = pd.DataFrame(columns=all_cells)
for chr in chr_list:
    chr_cnvs = h5_data["/cnvs/" + chr].value.transpose()
    chr_data = pd.DataFrame(chr_cnvs[:,:tot_cells], columns=all_cells)
    df_10ex = df_10ex.append(chr_data)


num_bin = df_10ex.count()
#############################################################
#                        Correlation                        #
#############################################################
corr = dict()
pvalue = dict()

for cell in cells:
    cellprofile_10ex = df_10ex.loc[:,(cell)].values
    cellprofile_ginkgo = df_ginkgo.loc[:,(cell)].values
    bad = np.in1d(cellprofile_10ex, range(0, -129), invert=True).reshape(cellprofile_10ex.shape)
    '''
    for pos, cn in enumerate(cellprofile_10ex):
        if cn < 0:
            if cn > -128:
                cellprofile_10ex[pos] = -cellprofile_10ex[pos]
    '''
    cellprofile_10ex = np.compress(bad, cellprofile_10ex)
    cellprofile_ginkgo = np.compress(bad, cellprofile_ginkgo)
    sp  = spearmanr(cellprofile_10ex, cellprofile_ginkgo)
    corr[cell] = sp[0]
    pvalue[cell] = sp[1]
    
#print('(Step 5/7) Plotting correlation distribution')

xaxis = []; yaxis = []
iaxis = []; jaxis = []
for k,v in corr.items():
    is_noisy = stats_df['is_noisy'].loc[stats_df['barcode'] == k].values[0]
    if is_noisy == 0:
        xaxis.append(k)
        yaxis.append(float(v))
    else:
        iaxis.append(k)
        jaxis.append(float(v))

plt.plot(range(len(yaxis)), yaxis, 'o', color='black', label='Non-noisy cells')
plt.xticks(range(len(xaxis)), list(xaxis))

plt.plot(range(len(jaxis)), jaxis, '*', color='red', label='Noisy cells')
plt.xticks(range(len(iaxis)), list(iaxis))

N = num_cells*0.01 #scale factor to avoid a huge png
plt.gca().margins(x=0)
plt.gcf().canvas.draw()
tl = plt.gca().get_xticklabels()
maxsize = max([t.get_window_extent().width for t in tl])
m = 0.2 # inch margin
s = maxsize/plt.gcf().dpi*N+2*m
margin = m/plt.gcf().get_size_inches()[0]
plt.gcf().subplots_adjust(left=margin, right=1.-margin, bottom=0.3)
plt.gcf().set_size_inches(s, plt.gcf().get_size_inches()[1])
plt.xticks([], [])
plt.gcf().suptitle('Spearman correlation - Sample: ' + sample + ', Cell type: ' + filter)
plt.gca().legend()

plt.savefig(outdir+"/spearman_noinf-int.png")
plt.close()

#print('(Step 6/7) Writing '+ 'pearson_corr.csv')
with open(outdir+'/spearman_noinf-int.csv', "w+") as f:
    f.write('cellid\tpearson corr coeff\tp value\t is noisy\n')
    for cell in cells:
        is_noisy = stats_df['is_noisy'].loc[stats_df["barcode"] == cell].values[0]
        if is_noisy == 0:
            f.write(cell+'\t'+ str(corr[cell])+'\t'+str(pvalue[cell])+'\t'+str(0)+'\n')
        else:
            f.write(cell+'\t'+ str(corr[cell])+'\t'+str(pvalue[cell])+'\t'+str(1)+'\n')

