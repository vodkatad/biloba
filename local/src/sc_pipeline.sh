#!/bin/bash

# run demultiplexer
snakemake -n demux.done

# run ginkgo analysis and draw plots
snakemake -n ginkgo.done

# run multi-sample analysis
snakemake -n multi_sample.done
