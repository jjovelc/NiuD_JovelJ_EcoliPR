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
### **4. Filtering variants**

- **Tool**: Biopet VCF Filter
- **Script**: `filter_vcf.sh`
- **Description**: Raw variants were filtered to retain high-quality calls using stringent quality thresholds.
- **Filtering Criteria**:
  - Minimum Sample Depth: ≥50 reads per sample
  - Minimum Quality Score: ≥30 (Phred-scaled)
  - Minimum Alternate Depth: ≥3 reads supporting the alternate allele

```bash
eval "$(conda shell.bash hook)"
conda activate call_variants_nano
for FILE in E.coli_NZ_HG941718.1.variants.vcf; do
  biopet-vcffilter -I ${FILE} \
                   -o ${FILE/.vcf/_filtered50-30-0.05.vcf} \
                   --minSampleDepth 50 \
                   --minQualScore 30 \
                   --minAlternateDepth 3
done
```

**Output**: `E.coli_NZ_HG941718.1.variants_filtered50-30-0.05.vcf`

---

### **5. Annotation of variants**

- **Tool**: snpEff
- **Script**: `run_snpEff.sh`
- **Description**: Filtered variants were functionally annotated to predict their effects on genes and proteins.
- **Database Build**: Custom snpEff database built from the reference genome and GFF annotation file.
- **GFF File Location**: The correct GFF file used by snpEff is located at `NZ_HG941718.1/genes.gff` within the snpEff data directory.

```bash
#!/usr/bin/bash

# Build custom snpEff database
java -jar ../snpEff.jar build -gff3 -v NZ_HG941718.1

# Annotate variants
java -jar ../snpEff.jar -v NZ_HG941718.1 \
     ../../E.coli_NZ_HG941718.1.variants_filtered50-30-0.05.vcf > \
     E.coli_NZ_HG941718.1.variants_filtered50-30-0.05_annotated.vcf
```

**Output**: `E.coli_NZ_HG941718.1.variants_filtered50-30-0.05_annotated.vcf`

---

### **6. Creation of genotype matrix**

- **Tool**: Python pandas
- **Script**: `create_matrix.py`
- **Description**: Annotated VCF file was converted into a numerical genotype matrix suitable for GWAS analysis. The script parses the VCF file and extracts genotype information for all samples across all variants.
- **Matrix Structure**:
  - Rows: SNP variants (identified as CHROM_POS)
  - Columns: Sample IDs (111 E. coli strains)
- **Genotype Encoding**:
  - 0/0 (homozygous reference) → 0
  - 0/1 (heterozygous) → 1
  - 1/1 (homozygous alternate) → 2
  - ./. (missing data) → -1

```python
python create_matrix.py E.coli_NZ_HG941718.1.variants_filtered50-30-0.05_annotated.vcf
```

**Output**: `genotype_matrix.csv` (tab-delimited matrix with variants as rows and 111 samples as columns)

### **7. Bayesian GWAS Analysis**

- **Tool**: [GWAS.BAYES (BICOSS)](https://cran.r-project.org/web/packages/GWAS.BAYES/)
- **Script**: `implement_bicoss.R`
- **Description**: A Bayesian genome-wide association study was performed for each phage resistance trait using BICOSS. The script processes the genotype matrix and phenotype data, calculates a kinship matrix, filters SNPs, and identifies statistically significant associations.

```r
# Load packages
library(GWAS.BAYES)

# Define input files
genotype_matrix <- "genotype_matrix.tsv"
traits_matrix <- "phage_traits.tsv"

# Run BICOSS with preprocessing, validation, and trait-wise GWAS
source("implement_bicoss.R")
```

---

### **8. Gene Name Mapping**

- **Tool**: R (`extract_gene_names.R`)
- **Script**: `extract_gene_names.R`
- **Description**: This step maps significant SNPs to gene names using the `genes.gff` file from the snpEff annotation. SNP positions are matched against coding sequences (CDS) in the GFF file, and nearby genes (within 2 kb) are assigned if no direct overlap is found.

```r
# Example usage
Rscript extract_gene_names.R
```
---

### **9. Visualization of SNP-Trait Associations**

- **Tool**: R (`ggplot2`, `dplyr`)
- **Script**: `plot_sign_SNPs.R`
- **Description**: This step generates a dot plot to visually summarize significant SNP-trait associations. Dot size reflects association strength (−log₁₀ p-value), and color indicates the direction and magnitude of effect size.

```r
# Example usage
Rscript plot_sign_SNPs.R
```

### 10. Generation of Variant Summary Table
- **Tool**: R
- **Script**: `generate_summary_table.R`
- **Description**: Creates a detailed summary table that combines SNPs, gene names, protein products, allele frequencies, and traits.

### Main results:

#### Dotplot including significant associations:

[Dotplot](data/second_iteration/SNP-trait-association-dotplot.pdf)

#### Summary table including all significant SNPs


