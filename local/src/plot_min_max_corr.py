import sys, os
import pandas as pd
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt

if len(sys.argv) < 6:
    print('Missing outpur: correct usage: ' + sys.argv[0] + 'pearson_correlation.csv segcopy_ginkgo cnvdata_10x per_cell_summary_metrics.csv outdir')
    exit(1)

if os.path.exists(sys.argv[1]) == False:
    print("FileNotFoundError: No such file or directory: '" + sys.argv[1] + "'")
    exit(1)

if os.path.exists(sys.argv[2]) == False:
    print("FileNotFoundError: No such file or directory: '" + sys.argv[2] + "'")
    exit(1)

if os.path.exists(sys.argv[3]) == False:
    print("FileNotFoundError: No such file or directory: '" + sys.argv[3] + "'")
    exit(1)

if os.path.exists(sys.argv[4]) == False:
    print("FileNotFoundError: No such file or directory: '" + sys.argv[4] + "'")
    exit(1)

if os.path.exists(sys.argv[5]) == False:
    print("FileNotFoundError: No such file or directory: '" + sys.argv[5] + "'")
    exit(1)  


correlation = sys.argv[1]
segcopy = sys.argv[2]
cnvdata = sys.argv[3]
cellmetrics = sys.argv[4]
outdir=sys.argv[5]

correlation_df = pd.read_csv(correlation, sep="\t")
segcopy_df = pd.read_csv(segcopy, sep="\t")
cellmetrics_df = pd.read_csv(cellmetrics, sep=",")

# find the max and the min correlation among cells
pearson_array = correlation_df.loc[:,('pearson corr coeff')]

max = pearson_array.max()
min = pearson_array.min()

# find the position of the min and the max correlation in the dataframe
#idxmax = np.where(pearson_array == max)[0]
#idxmin = np.where(pearson_array == min)[0]

# retrieve the barcodes of the cells with max and min correlation
barcodemax=correlation_df.loc[correlation_df['pearson corr coeff'] == max]['cellid'].values[0]
barcodemin=correlation_df.loc[correlation_df['pearson corr coeff'] == min]['cellid'].values[0]

# retrieve ginkgo cnv profiles for the selected cells
cnv_max_ginkgo = segcopy_df.loc[:,(barcodemax)].values
cnv_min_ginkgo = segcopy_df.loc[:,(barcodemin)].values

idmax = cellmetrics_df.loc[cellmetrics_df['barcode'] == barcodemax + '-1']['cell_id'].values[0]
idmin = cellmetrics_df.loc[cellmetrics_df['barcode'] == barcodemin + '-1']['cell_id'].values[0]

# retrieve 10x cnv profiles for the selected cells

h5_f = h5.File(sys.argv[3])
chr_list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
            "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
cnv_max_10ex = []
cnv_min_10ex = []
for chr in chr_list:
    chr_cnvs = h5_f["/cnvs/" + chr].value
    chr_data_max = chr_cnvs.transpose()[:,idmax]
    chr_data_min = chr_cnvs.transpose()[:,idmin]
    cnv_max_10ex = np.concatenate([cnv_max_10ex, chr_data_max])
    cnv_min_10ex = np.concatenate([cnv_min_10ex, chr_data_max])

# change the sign of 10x negative cnvs (inferred) and mask cnvs=-128 
bad = np.in1d(cnv_max_10ex, np.arange(-1,-129, -1), invert=True).reshape(cnv_max_10ex.shape)
#for idx, v in enumerate(cnv_max_10ex):
    #if v < 0:
        #if v > -128:
            #cnv_max_10ex[idx] = -cnv_max_10ex[idx]
cnv_max_10ex = np.compress(bad, cnv_max_10ex)
cnv_max_ginkgo = np.compress(bad, cnv_max_ginkgo)
            
bad = np.in1d(cnv_min_10ex, np.arange(-1, -129, -1), invert=True).reshape(cnv_min_10ex.shape)
for idx, v in enumerate(cnv_min_10ex):
    if v < 0:
        if v > -128:
            cnv_min_10ex[idx] = -cnv_min_10ex[idx]
cnv_min_10ex = np.compress(bad, cnv_min_10ex)
cnv_min_ginkgo = np.compress(bad, cnv_min_ginkgo)

#print("len(cnv_max_10ex): " + str(len(cnv_max_10ex)))
#print("len(cnv_max_ginkgo): " + str(len(cnv_max_10ex)))


# plot cnv profiles of the cell with max correlation
plt.hexbin(cnv_max_ginkgo, cnv_max_10ex, gridsize=(1500, 1500), cmap='inferno' )
plt.xlabel("Ginkgo cnvs")
plt.ylabel("10x cnvs")
plt.gcf().suptitle('Ginkgo CNVS vs 10x CNVS - Max Pearson Coeff', fontsize=14)
#plt.xlim(2002, 2008)
#plt.ylim(0, 4500)
#props = dict(boxstyle='round', facecolor='white', alpha=0.5)
#s = "cellid: " + barcodemax + "\nPearson Corr Coeff: " + str(max)
#leg = plt.figlegend([a, b, c, d], ['1','2', '3', '4'], loc=(0.85, 0.65))
#plt.figtext(0.5, 0, s, wrap=True, horizontalalignment='center', fontsize=10)
#plt.subplots_adjust(left=0.25)
plt.savefig(outdir + '/max_pearson_scatter_hexbin.png')
plt.clf()

plt.hexbin(cnv_min_ginkgo, cnv_min_10ex, gridsize=(1500, 1500), cmap='inferno')
plt.xlabel("Ginkgo cnvs")
plt.ylabel("10x cnvs")
plt.gcf().suptitle('Ginkgo CNVS vs 10x CNVS - Min Pearson Coeff', fontsize=14)
#props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#s = "cellid: " + barcodemin + "\nPearson Corr Coeff: " + str(min)
#plt.text(0.05, 0.95, s, transform=plt.gca().transAxes, fontsize=10,
#        verticalalignment='top', bbox=props)

plt.savefig(outdir+'/min_pearson_scatter_hexbin.png')
plt.clf()



