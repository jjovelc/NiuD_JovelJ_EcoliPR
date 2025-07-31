# SNP Discovery and GWAS Pipeline for Phage Resistance in *E. coli*

(Analyses were revisited for a grant application).

This repository documents the variant calling and genome-wide association study (GWAS) pipeline used to identify significant non-synonymous SNPs (nsSNPs) associated with phage resistance in 111 cultured *Escherichia coli* isolates.

## Overview

A total of 111 *E. coli* strains were sequenced, aligned to a reference genome, and analyzed to discover genetic variants. Variants were filtered, annotated, and subjected to association analysis against a phage resistance phenotype using Bayesian GWAS (BICOSS). This document details the computational steps used from raw sequencing data to the identification of statistically significant SNPs.

---

## Pipeline Steps

### 1. Input Data

- **Samples**: 111 cultured *E. coli* strains
- **Reference Genome**: `NZ_HG941718.1` (used for both alignment and variant calling)
- **Annotation File**: GFF file generated with Prokka from `NZ_HG941718.1`
- **phage_traits.tsv**: A matrix with one E.coli strain per row and a Phage per column. Data is a matrix with 0 or 1 values depending if a bacteria was or not susceptible to each phage.

---

### 2. Read Alignment

- **Tool**: [BWA](http://bio-bwa.sourceforge.net/)
- **Script**: `align_against_assemblies.sh`
- **Description**: Paired-end reads were aligned to the reference genome.
- **Output**: Sorted BAM files

---

### 3. Variant Calling

- **Tool**: [FreeBayes](https://github.com/freebayes/freebayes)
- **Script**: `run_freebayes.sh`
- **Description**: Multi-sample variant calling was performed on the BAM files.
- **Output**: `variants.vcf`

---

### 4. Filtering variants




