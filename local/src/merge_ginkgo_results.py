import pandas as pd
import argparse
import json
# - s results.txt,SegFixed,SegBreaks s results.txt,SegFixed,SegBreaks [-s3...]

#output:
# - results_merged.txt
# - SegFixed_merged
# - SegBreaks_merged

def required_length(nmin):
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not nmin<=len(values):
                msg='argument "{f}" requires at least {nmin} arguments'.format(
                        f=self.dest,nmin=nmin)
                raise argparse.ArgumentTypeError(msg)
            
            setattr(args, self.dest, values)
    return RequiredLength


parser = argparse.ArgumentParser(description="Ginkgo analysis results on multiple samples merging")
parser.add_argument("sample", metavar='sample_name=Results.txt,SegFixed,SegBreaks', action=required_length(2),
        help='Sample name followed by the path to the result files produced by Ginkgo in the order indicated. At least two samples must be provided.',
        nargs='+', type=str)
parser.add_argument("outdir", metavar='path/to/out/dir', action='store',
        help='Path to the desired output directory where the merged files have to be stored',
        nargs=1, type=str)
args=parser.parse_args()

samples = args.sample
outdir = args.outdir[0]

samples_dict = {}
for sample in samples:
    sample_name = sample.split("=")[0]
    files = sample.split("=")[1].split(",")
    
    samples_dict[sample_name] = {
            'results' : files[0],
            'segfixed' : files[1],
            'segbreaks' : files[2],
            }


'''
    Create results_merged.txt file and rename cells
'''

sample_cells_dict = {}
merged_results = outdir+"/results_merged.txt"
with open(merged_results, "w+") as out_f:
    cellcounter = 1
    header = "Sample\tSoSPredictedPloidy(Top5)\tErrorInSoSApproach(Top5)\tCopyNumber(SoS)"
    out_f.write(header+'\n')
    for sample_name, files in samples_dict.items():
        sample_cells_dict[sample_name] = {}
        with open(files['results'], 'r') as in_f:
            isheader = True
            for line in in_f:
                if isheader == False: #skip header
                    '''
                        === 1 ===
                    '''
                    newline = '=== ' + str(cellcounter) + '==='
                    #out_f.write(newline+'\n')
                    
                    '''
                        [1] ...
                    
                    out_f.write(next(in_f))
                    '''

                    '''
                        SAMPLE ....
                    '''
                    fields = next(in_f).split('\t')
                    cellid= fields[0]
                    sample_cells_dict[sample_name][cellid] = str(cellcounter - 1)
                    newline = str(cellcounter - 1) + '\t' + fields[1] + '\t' + fields[2] + '\t' + fields[3]
                    out_f.write(newline)
                    
                    cellcounter+=1
                else:
                    #the header line resets the flag
                    isheader = False

'''
    Create samples.json file, a dictionary storing the
    mapping among sample_names-cellids-new_cellids
'''

samples_json = outdir + "/samples.json"
with open(samples_json, 'w+') as out_f:
    json.dump(sample_cells_dict, out_f)

'''
    Create SegFixed_merged
'''
segfixed_merged = outdir + '/SegFixed_merged'

columns_copied = True #to force to not c
df_merged = pd.DataFrame()
for sample_name, files in samples_dict.items():
    df = pd.read_csv(files['segfixed'], sep="\t")
    #print(df)
    if columns_copied == False:
        df_merged['CHR'] = df['CHR']
        df_merged['START'] = df['START']
        df_merged['END'] = df['END']
        columns_copied = True
    for column in df:
        if column != 'CHR':
            if column != 'START':
                if column != 'END':
                    if column in sample_cells_dict[sample_name].keys(): #to avoid unnamed columns
                        cellid = sample_cells_dict[sample_name][column] #retrieve the new cellid
                        df_merged[cellid] = df[column]

df_merged.to_csv(segfixed_merged, sep="\t", index=False)

'''
    Create SegBreaks_merged
'''
segfixed_merged = outdir + '/SegBreaks_merged'

columns_copied = True
df_merged = pd.DataFrame()
for sample_name, files in samples_dict.items():
    df = pd.read_csv(files['segbreaks'], sep="\t")
    #print(df)
    if columns_copied == False:
        df_merged['CHR'] = df['CHR']
        df_merged['START'] = df['START']
        df_merged['END'] = df['END']
        columns_copied = True
    for column in df:
        if column != 'CHR':
            if column != 'START':
                if column != 'END':
                    if column in sample_cells_dict[sample_name].keys(): #to avoid unnamed columns
                        cellid = sample_cells_dict[sample_name][column] #retrieve the new cellid
                        df_merged[cellid] = df[column]

df_merged.to_csv(segfixed_merged, sep="\t", index=False)















