import sys, os
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description="Removes normal cells and recomputes phylogenetic trees")

parser.add_argument("ploidy", metavar='low_thresh high_thresh', action='store', 
        help='Lower and higher limits to the ploidy range to be discarded.', 
        nargs=2, type=float)
parser.add_argument("results", metavar='path/to/results.txt', action='store',
        help='Ginkgo results.txt file',
        nargs=1, type=str)
parser.add_argument("segcopy", metavar='path/to/SegCopy', action='store',
        help='Ginkgo SegCopy file.',
        nargs=1, type=str)
parser.add_argument("outdir", metavar='path/to/out_directory', action='store',
        help='Path to output directory.',
        nargs=1, type=str)

args=parser.parse_args()
low_thresh = args.ploidy[0]
high_thresh = args.ploidy[1]
results = args.results[0]
segcopy = args.segcopy[0]
outdir = args.outdir[0]

ok_cells_f = outdir + "/ok_cells.txt"
ko_cells_f = outdir + "/ko_cells.txt"

'''
    Filter out rows (cells) where CopyNumber is equal to 2.

'''

# Sample    CopyNumber(SoS) SoSPredictedPloidy(Top5)          ErrorInSoSApproach(Top5)

df = pd.read_csv(results, sep="\t")
df = df.dropna()

df_f = df.loc[(df["Copy_Number"] < low_thresh) | (df["Copy_Number"] > high_thresh)]
# filtered_results.txt -> results in the desired threshold
df_f.to_csv(outdir+"/filtered_results.txt", sep="\t", index=False)
filtered_cells = df_f.Sample.values
columns = np.append(filtered_cells, ['CHR', 'START', 'END'])

segcopy_df = pd.read_csv(segcopy, sep="\t", usecols = lambda column : column in columns)
segcopy_df.to_csv(outdir+"/filtered_SegCopy", sep="\t", index=False)
with open(ok_cells_f, "w+") as f:
    for cell in filtered_cells:
        f.write(cell+'\n')




