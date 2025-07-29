#!/usr/bin/bash

INFILE=$1
OUTFILE=$2

java -jar snpEff.jar ann \
    -no-downstream \
    -no-upstream \
    -no-intergenic \
    NZ_HG941718.1 \
    "$INFILE" \
    > "$OUTFILE"
