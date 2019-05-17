import sys
import os  
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

if len(sys.argv) < 3:
    print('Missing arguments --- Correct usage: ' + argv[0] + 'SegCopy out_path')
    exit(1)

if os.path.isfile(sys.argv[1]) == False:
    print('FileNotFoundError: No such file or directory:' + sys.argv[1])
    exit(1)

if os.path.isdir(sys.argv[2]) == False:
    print('FileNotFoundError: No such file or directory:' + sys.argv[2])
    exit(1)

df = pd.read_csv(sys.argv[1], sep="\t")
chromosomes = df["CHR"].unique()

for chrom in chromosomes:
    chr_df = df.loc[df['CHR'] == chrom]
    chr_df.index = chr_df["START"]
    chr_df = chr_df.drop(["CHR", "START", "END"], axis=1)
    plt.figure(figsize=(18.53, 9.91), dpi=100)
    sns.heatmap(chr_df.T, vmax=6)
    print('Creating ' + chrom + '_heatCN.png')
    plt.legend(loc='best')
    plt.savefig(sys.argv[2] + '/' + chrom + '_heatCN.png')
    plt.clf()
