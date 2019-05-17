import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Density plot of cnv events associated to a gene.")
parser.add_argument("file", metavar='genes_cnvs', action='store',
                help='File containing the cnvs inferred on genes.',
                        nargs=1, type=str)
parser.add_argument("gene", metavar='gene', action='store',
                        help='Gene to be plotted.', nargs=1, type=str)

parser.add_argument("outfile", metavar='path/to/out.png', action='store',
                help='Output file.',
                        nargs=1, type=str)
args=parser.parse_args()

genes_cnvs = args.file[0]
gene = args.gene[0] 
outfile = args.outfile[0]

df = pd.read_csv(genes_cnvs, sep='\t', header=None)
cnvs = df[df[3] == gene][4]
cnvs = np.asarray(cnvs[cnvs.index[0]].split(';')).astype(int)

sns.set_style('whitegrid')
sns.kdeplot(cnvs, bw=0.5)
plt.savefig(outfile)
