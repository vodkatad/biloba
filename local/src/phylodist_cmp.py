import sys
import json
import numpy as np
import pandas as pd
import dendropy
import argparse
from dendropy import treecalc
import pandas as pd


parser = argparse.ArgumentParser(description="Multi sample phylogenetic distance")
parser.add_argument("newick", metavar='multi-sample-tree.newick', action='store',
                help='File containing the aggregate phylogenetic tree in newick format',
                        nargs=1, type=str)
parser.add_argument("json", metavar="samples.json", action='store',
        help='Samples info json file', nargs=1, type=str)
parser.add_argument("outdir", metavar='path/to/out/dir', action='store',
                help='Path to the desired output directory where the merged files have to be stored',
                        nargs=1, type=str)
args=parser.parse_args()

clust_f = args.newick[0]
samples_f = args.json[0]
outdir = args.outdir[0]


def phylogenetic_distance_matrix_from_file(file_path):
    t = dendropy.Tree.get_from_path(file_path, "newick")
       
    #patristic distance matrix
    pdm = t.phylogenetic_distance_matrix().as_data_table()
    dm = pd.DataFrame(index=pdm._row_names, columns=pdm._column_names)
    for row in pdm._row_names:
        dm.loc[row] = pdm._data[row]    
    return dm

def split_distance_matrix(dm, row_list, column_list):
    dm_ = pd.DataFrame()
    dm_ = dm[column_list]
    dm_ = dm_.ix[row_list]
    return dm_

def intra_sample_mean_pairwice_distance(dm):
    distances = []
    for i in range(0, len(dm.values)):
        for j in range(i+1, len(dm.values)):
            distances.append(dm.values[i,j])
    return sum(distances)/len(distances)


def inter_sample_mean_pairwice_distance(dm):
    distances = []
    for row_name in dm.index:
        row = dm.loc[row_name]
        for col_name in dm.columns:
            distances.append(row[col_name])
    return sum(distances)/len(distances)

#print("Generating distance matrix")
dm = phylogenetic_distance_matrix_from_file(clust_f)
    
#read json file
with open(samples_f) as f:  
    samples_info = json.load(f)

#print("Generating single sample distance matrices")
#generate one distance matrix for each sample
samples = list(samples_info.keys())
sample_dms = {}
for sample in samples:
    id_list = list(samples_info[sample].values())
    s_dm = split_distance_matrix(dm, id_list, id_list)
    sample_dms[sample] = s_dm

#print("Computing intra-sample mean pairwise distance")
# compute intra-sample mean pairwise distance
sample_mpd = pd.DataFrame(columns=['sample','intra-th'])
for sample in samples:
    s_dm = sample_dms[sample]
    sample_mpd = sample_mpd.append({'sample':sample, 'intra-th':intra_sample_mean_pairwice_distance(s_dm)}, ignore_index=True)
    #print(sample + " -- mean pairwise distance: " + str(sample_mpd[sample]))
sample_mpd.to_csv(outdir+"/intra-sample-het.csv", sep="\t", index=False)

#print("Generating sample-couples distance matrices")
#generate one distance matrix for each couple of samples
#samples = samples_info['sample_names']
couple_dms = {}
for sample in samples:
    couple_dms[sample] = {}
for i in range(0, len(samples)-1):
    row_names = list(samples_info[samples[i]].values())
    for j in range(i+1, len(samples)):
        col_names = list(samples_info[samples[j]].values())
        c_dm = split_distance_matrix(dm, row_names, col_names)
        couple_dms[samples[i]][samples[j]] = c_dm
        couple_dms[samples[j]][samples[i]] = c_dm
    
#print("Computing sample-couples mean pairwise distance")
# compute inter-sample mean pairwise distance
couple_mpd = pd.DataFrame(columns=['sample1','sample2','inter-th'])
for i in range(0, len(samples)-1):
    for j in range(i+1, len(samples)):
        c_dm = couple_dms[samples[i]][samples[j]]
        couple_mpd = couple_mpd.append({'sample1':samples[i], 'sample2':samples[j], 'inter-th':inter_sample_mean_pairwice_distance(c_dm)}, ignore_index=True)
        #couple_mpd[samples[i]][samples[j]] = inter_sample_mean_pairwice_distance(c_dm)
        #couple_mpd[samples[j]][samples[i]] = inter_sample_mean_pairwice_distance(c_dm)
        #print(samples[i] + " - " + samples[j] +  " -- mean pairwise distance: " + str(couple_mpd[samples[i]][samples[j]]))
couple_mpd.to_csv(outdir+"/inter-samples-het.csv", sep="\t", index=False)
