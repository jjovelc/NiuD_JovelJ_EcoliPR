# 	GWAS analysis to identify E. coli genes or SNPS associated with phage resistance

# SNP Discovery and GWAS Pipeline for Phage Resistance in *E. coli*

This repository documents the variant calling and genome-wide association study (GWAS) pipeline used to identify significant non-synonymous SNPs (nsSNPs) associated with phage resistance in 111 cultured *Escherichia coli* isolates.

## Overview

A total of 111 *E. coli* strains were sequenced, aligned to a reference genome, and analyzed to discover genetic variants. Variants were filtered, annotated, and subjected to association analysis against a phage resistance phenotype using Bayesian GWAS (BICOSS). This document details the computational steps used from raw sequencing data to the identification of statistically significant SNPs.

---

## Pipeline Steps

### 1. Input Data

- **Samples**: 111 cultured *E. coli* strains  
- **Reference Genome**: `NZ_HG941718.1` (used for both alignment and variant calling)  
- **Annotation File**: GFF file generated with Prokka from `NZ_HG941718.1`  

---

### 2. Read Alignment

- **Tool**: [BWA](http://bio-bwa.sourceforge.net/)  
- **Description**: Paired-end reads were aligned to the reference genome.  
- **Output**: Sorted BAM files  

---

### 3. Variant Calling

- **Tool**: [FreeBayes](https://github.com/freebayes/freebayes)  
- **Description**: Multi-sample variant calling was performed on the BAM files.  
- **Output**: `variants.vcf`  

---

### 4. Variant Filtering

- **Tool**: [biopet-vcffilter](https://github.com/biopet/biopet)  
- **Parameters**:
  - Minimum read depth: `50`
  - Minimum quality score: `>= 30`
  - Minimum alternate allele depth: `5%`  
- **Output**: `filtered.vcf`  

---

### 5. Variant Annotation

- **Tool**: [snpEff](https://pcingola.github.io/SnpEff/)  
- **Description**: SNPs were annotated using a snpEff database generated with Prokka annotation.  
- **Input**: `filtered.vcf`, Prokka-derived GFF  
- **Output**: `annotated.vcf`  

---

### 6. Extraction of Non-Synonymous SNPs

- **Tool**: [VCFtools](https://vcftools.github.io/)  
- **Description**: Non-synonymous SNPs (nsSNPs) were filtered out of the annotated VCF.  
- **Output**: A list of nsSNPs used in the GWAS  

---

### 7. GWAS Analysis

- **Tool**: [BICOSS](https://github.com/maljovec/bicoss) (Bayesian GWAS implemented in GAWS.BAYES)  
- **Phenotype**: Phage resistance  
- **Input**: Matrix of nsSNPs vs samples, with binary resistance phenotype  
- **Statistical Threshold**: FDR = 0.05  
- **Output**: Significant SNPs associated with phage resistance  

---

## Output Files

- `BAM files`: Aligned sequencing reads  
- `variants.vcf`: Raw SNPs called across samples  
- `filtered.vcf`: Quality-filtered variants  
- `annotated.vcf`: SNPs annotated with effect information  
- `nsSNPs list`: Subset of SNPs with non-synonymous effects  
- `significant SNPs`: SNPs significantly associated with phage resistance  

---

## Dependencies

Please ensure the following software is installed and available in your environment:

- BWA  
- FreeBayes  
- biopet-vcffilter  
- snpEff  
- VCFtools  
- BICOSS (GAWS.BAYES)  

---

## Notes

- The reference genome and annotation must be consistent across all steps.  
- Filtering thresholds (e.g., depth, QScore) can be adjusted as needed depending on your sequencing depth and sample quality.  
- The pipeline assumes samples are phenotyped for resistance/sensitivity to a given phage.  


The analytical pipeline is described in [Fig. 1]("data/Fig.1_workflow.png")

Figure 1. Workflow used for identification of SNPs in E. coli genomic libraries and GWAS analysis.

