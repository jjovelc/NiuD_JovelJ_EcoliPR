#!/usr/bin/bash

INFILE=$1
OUTFILE=$2

bcftools query -f '[%CHROM:%POS\t%SAMPLE\t%GT\n]' "$INFILE" > "$OUTFILE"
