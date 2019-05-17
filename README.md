# biloba

A pipeline which combines the advantage of transparently dealing with 10x sequencing data of Cell Ranger DNA, 
with the openness and flexibility of Ginkgo, used in a stand-alone fashion, in order to perform a _multi-sample_
single-cell CNV analysis, on _large-scale_ datasets

## Installation

Requirements (tested versions in parentheses):

- Snakemake (5.4.4)
- python (3.5.3)
- gcc (5.4.0)
- R (3.5.1)

Python libraries:
- numpy
- pandas
- scipy
- seaborn
- argparse
- matplotlib

After having cloned the repositories and installed the requirements cellranger-dna and its
appropriate reference genomes tarballs need to be downloaded from the 10x website 
(https://support.10xgenomics.com/single-cell-dna/software/downloads/latest).
A bash script to perform the whole analysis on public datasets is available:
`dataset/publicdataset/cellrangerdna/sc_pipeline`
In this same directory users should put the results from cellranger-dna from public datasets.
For your own 10x experiments the pipeline is in `dataset/cellrangerdna`, the directory
where the fastq files can be found has to be set up in the `conf.sk` file.
