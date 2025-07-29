#!/usr/bin/bash

INFILE=$1
OUTFILE=$2

bcftools view -i 'INFO/ANN[0] !~ "synonymous_variant"' "$INFILE" > "$OUTFILE"
