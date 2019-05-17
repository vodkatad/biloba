import sys, os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

parser = argparse.ArgumentParser(description="Cells mean ploidy density plot (ginkgo)")
parser.add_argument("results", metavar='path/to/results.txt', action='store',
                help='Ginkgo results.txt file',
                        nargs=1, type=str)
parser.add_argument("tool", metavar='ginkgo/10x', choices=['ginkgo', '10x'], action='store',
                        help='Tool used to produce the data to plot',
                                                nargs=1, type=str)

parser.add_argument("sample", metavar='Sample1', action='store',
        help='Sample name',
        nargs=1, type=str)


parser.add_argument("outdir", metavar='path/to/out/dir', action='store',
                help='Path to the output directory',
                        nargs=1, type=str)
args=parser.parse_args()

results = args.results[0]
tool = args.tool[0]
outdir = args.outdir[0]

if tool == 'ginkgo':
    df = pd.read_csv(results, sep="\t")
    df = df.dropna()
    keyword = 'Copy_Number'
else:
    df = pd.read_csv(results)
    keyword = 'mean_ploidy'

min_cn = np.min(df[keyword].values)
max_cn = np.max(df[keyword].values)

df[keyword].plot.kde(bw_method=0.05)

plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.xticks(np.arange(min_cn, max_cn+1, 0.1))
plt.xlim(min_cn, max_cn)

plt.grid()

plt.gcf().set_size_inches(37,21)
plt.gcf().subtitle("Mean CN density plot - Tool = " + tool + ", Sample = " + sample)
plt.savefig(outdir+"/density_CN.png", dpi=100)
