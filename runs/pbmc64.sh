#!/bin/bash
mkdir -p output/
CMD="julia -p 64 scripts/runPIDCParallel.jl -e . -r 4 --hdf5 ./data/pbmc/adata.20K.5K.h5ad ./output/pidc_64p.csv"
echo $CMD
$CMD 2>&1 | tee output/pidc.log
$CMD 2>&1 | tee -a output/pidc.log
$CMD 2>&1 | tee -a output/pidc.log
