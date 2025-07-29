#!/usr/bin/bash

eval "$(conda shell.bash hook)"
conda activate prokka

INFILE=$1

prokka --centre X --compliant --cpus 8 --outdir "${INFILE/.fasta/}"_prokka  --prefix "${INFILE/.fasta/}" $INFILE
