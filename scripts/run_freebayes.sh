#!/usr/bin/bash

# Run FreeBayes in parallel
# The genome is divided into 100,000-base regions for parallel processing.
# Using 48 threads for high-performance processing.
freebayes-parallel \
    <(fasta_generate_regions.py \
        /work/vetmed_data/jj/projects/dongyanNiu/ST131_Data_2024/all_fastq/NZ_HG941718.1.fasta.fai 100000) \
    48 \
    -f /work/vetmed_data/jj/projects/dongyanNiu/ST131_Data_2024/all_fastq/NZ_HG941718.1.fasta \
    --ploidy 1 \
    --min-alternate-total 2 \
    /work/vetmed_data/jj/projects/dongyanNiu/ST131_Data_2024/all_fastq/bamfiles/*.bam \
    > /work/vetmed_data/jj/projects/dongyanNiu/ST131_Data_2024/all_fastq/bamfiles/E.coli_NZ_HG941718.1.variants.vcf
