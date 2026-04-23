#!/bin/bash

mkdir -p data/

wget -O data/pbmc20K.zip https://zenodo.org/records/19670636/files/pbmc20K.zip?download=1

unzip data/pbmc20K.zip -d data/
