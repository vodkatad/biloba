import pandas as pd
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster import hierarchy
import matplotlib.colors as colors
import json

def required_length(nmin):
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin<=len(values):
                msg='argument "{f}" requires at least {nmin} arguments'.format(
                        f=self.dest,nmin=nmin)
                raise argparse.ArgumentTypeError(msg)
            
            setattr(args, self.dest, values)
    return RequiredLength

def linkage_to_newick(Z, labels):
    """
    Input :  Z = linkage matrix, labels = leaf labels
    Output:  Newick formatted tree string
    """
    tree = hierarchy.to_tree(Z, False)
    def buildNewick(node, newick, parentdist, leaf_names):
        if node.is_leaf():
            return "{leaf_name}:{dist}{newick}".format(newick=newick, leaf_name=leaf_names[node.id], dist=(parentdist - node.dist)/2)
        else:
            if len(newick) > 0:
                newick = "):{dist}{newick}".format(newick=newick, dist=(parentdist - node.dist)/2)
            else:
                newick = ");"
            newick = buildNewick(node.get_left(), newick, node.dist, leaf_names)
            newick = buildNewick(node.get_right(), ",{newick}".format(newick=newick), node.dist, leaf_names)
            newick = "({newick}".format(newick=newick) 
            return newick
    return buildNewick(tree, "", tree.dist, labels)

parser = argparse.ArgumentParser(description="Ginkgo analysis results on multiple samples merging")
parser.add_argument("sample", metavar='sample_name=SegCopy', action=required_length(2),
        help='Sample name followed by the path to the CN file produced by Ginkgo. At least two samples must be provided.',
        nargs='+', type=str)
parser.add_argument("outdir", metavar='path/to/out/dir', action='store',
        help='Path to the desired output directory where the merged files have to be stored',
        nargs=1, type=str)
args=parser.parse_args()

samples = args.sample
outdir = args.outdir[0]

samples_dict = {}
for sample in samples:

    fields = sample.split("=")
    sample_name = fields[0]
    
    samples_dict[sample_name] = {
            'segcopy' : fields[1]
            }
samples = samples_dict.keys()

'''
    Create SegCopy_merged
'''
segcopy_merged = outdir + '/SegCopy_merged'
sample_cells_dict = {}
cn_df = pd.DataFrame()
cellcounter = 1
for sample_name, files in samples_dict.items():
    df = pd.read_csv(files['segcopy'], sep="\t",  usecols = lambda column : column not in ['CHR','START', 'END', 'Unnamed: 103', 'Unnamed: 113'])
    sample_cells_dict[sample_name] = {}
    for column in df:
        cellid = str(cellcounter - 1)
        sample_cells_dict[sample_name][column] = cellid
        cn_df[cellid] = df[column]
        cellcounter += 1
'''
    Create samples.json file, a dictionary storing the
    mapping among sample_names-cellids-new_cellids
'''

samples_json = outdir + "/samples.json"
with open(samples_json, 'w+') as out_f:
    json.dump(sample_cells_dict, out_f)


'''
    Generate multi-sample dendrogram and heatmap
'''
sns.set(color_codes=True)
cn_df = cn_df.transpose()

#add sample column
cn_df.index.name = 'cell'
cn_df['sample'] = ''

'''
iterate on rows, and check in which sample.values() ends up the index value (cellid) 
'''
for idx, row in cn_df.iterrows():
    for sample in samples:
        if idx in sample_cells_dict[sample].values():
            cn_df.loc[idx, 'sample'] = sample
#multilevel_index
cn_df.set_index([cn_df.index, 'sample'], inplace=True)

#http://dawnmy.github.io/2016/10/24/Plot-heatmaap-with-side-color-indicating-the-class-of-variables/

# Create a categorical palette to identify the samples
sample_pal = sns.color_palette("bright", 2)
#flatui = ["#ff3333", "#112be7"]
#sample_pal = sns.color_palette(flatui)
sample_lut = dict(zip(map(str, samples), sample_pal))
# Convert the palette to vectors that will be drawn on the side of the matrix
allsamples = cn_df.index.get_level_values('sample')

cn_df = cn_df.reset_index(level=1, drop=True)
sample_colors=pd.Series(allsamples, index=cn_df.index).map(sample_lut)

#norm=MidpointNormalize(vmin=0., midpoint=2., vmax=12.)

cbar_kws={"ticks":np.arange(0,6,1)}
h = sns.clustermap(cn_df, method='complete', metric='euclidean', 
        col_cluster=False, 
        yticklabels = False,
        row_colors=sample_colors, 
        cmap='RdBu_r',
        vmin=0, vmax=5,
        center=2,
        #norm=norm,
        cbar_kws=cbar_kws)

plt.gcf().suptitle('Multi-sample heatmap - Samples = ' + ', '.join(str(x) for x in samples)) 
plt.gcf().set_size_inches(37, 21)
plt.savefig(outdir+'/multi_sample_heatmap.png')                                  


'''
    Save multi-sample newick tree
'''

linkage = h.dendrogram_row.linkage
labels = h.dendrogram_row.dendrogram['leaves']
newick = linkage_to_newick(linkage,labels)

with open(outdir+'/multi_sample_tree.newick', 'w+') as f:
        f.write(newick)

