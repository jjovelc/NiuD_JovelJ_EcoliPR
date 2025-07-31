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

| Gene Name | Protein Product | Protein ID | Reference Genome Mutation Site | Reference Allele | SNP Allele | Percentage of isolates carrying this SNPs | Traits Associated |
|-----------|-----------------|------------|--------------------------------|------------------|------------|-------------------------------------------|-------------------|
| JKKDPAEF_00206 | IS66 family transposase ISEc23 | JKKDPAEF_00206 | 236104 | A | C | 1.8% | p186B |
| JKKDPAEF_00297 | hypothetical protein | JKKDPAEF_00297 | 318068 | GG | GA | 1.9% | p630P |
| JKKDPAEF_00890 | hypothetical protein | JKKDPAEF_00890 | 959719 | ACGCACA | GCGTAGC | 5.9% | p261P |
| agp | Glucose-1-phosphatase | JKKDPAEF_01073 | 1129589 | CCT | ACC | 29.8% | p331F |
| JKKDPAEF_01225 | hypothetical protein | JKKDPAEF_01225 | 1275120 | GTTATCTGCGTCA | GTTATCTGCATCA | 9.7% | p630P |
| JKKDPAEF_01379 | hypothetical protein | JKKDPAEF_01379 | 1418185 | GT | GC | 22.8% | p210B |
| JKKDPAEF_01379 | hypothetical protein | JKKDPAEF_01379 | 1418533 | A | T | 16.7% | p331F |
| JKKDPAEF_01400 | hypothetical protein | JKKDPAEF_01400 | 1440661 | GG | GA | 4.5% | p112PL |
| JKKDPAEF_01407 | hypothetical protein | JKKDPAEF_01407 | 1443597 | AAATTTTCCGTTGGCTGGAT | AAATTTTCCGCTGGGTAGAG | 3.6% | p069PL |
| JKKDPAEF_01456 | hypothetical protein | JKKDPAEF_01456 | 1486800 | AAGTTTTTCACTTCC | AAGTTATTTACTTCC | 0.9% | p630P |
| JKKDPAEF_01457 | hypothetical protein | JKKDPAEF_01457 | 1487649 | TTTCCCTAAGGCA | TTTTCCTGAAGCA | 0.9% | p186F |
| ompX_3 | Outer membrane protein X | JKKDPAEF_02111 | 2143638 | GCCCTTCAGGTCATCACTG | GCCCTTCAGGTCATCGCTG | 6% | p630P |
| JKKDPAEF_02112 | hypothetical protein | JKKDPAEF_02112 | 2144626 | ATGTTTG | ATGCTTA | 17.6% | p572FHCM |
| JKKDPAEF_02112 | hypothetical protein | JKKDPAEF_02112 | 2144902 | ATAA | ATGG | 26% | p572PLE |
| flk | Flagellar regulator flk | JKKDPAEF_02552 | 2630326 | G | A | 16.7% | p112PL |
| JKKDPAEF_03185 | IS3 family transposase IS2 | JKKDPAEF_03185 | 3324323 | G | A | 7.5% | p331F |
| JKKDPAEF_03197 | hypothetical protein | JKKDPAEF_03197 | 3335285 | GCCT | GCCA | 7.3% | p186B |
| JKKDPAEF_03219 | hypothetical protein | JKKDPAEF_03219 | 3358094 | TATCATCA | CATAACCA | 2.7% | p186P |
| kpsD | Polysialic acid transport protein KpsD | JKKDPAEF_03227 | 3367400 | AGGAGTTG | AGGGGTTG | 82.2% | p630P |
| JKKDPAEF_03229 | hypothetical protein | JKKDPAEF_03229 | 3369574 | AAAAGTT | AAAAATC | 78.5% | p069PL |
| JKKDPAEF_03229 | hypothetical protein | JKKDPAEF_03229 | 3370025 | ACACG | ACACT | 50.5% | p526RVM, p186B |
| JKKDPAEF_03230 | hypothetical protein | JKKDPAEF_03230 | 3370198 | GCAA | ACAA | 50.5% | p266F, p204B, p069FHC |
| JKKDPAEF_03230 | hypothetical protein | JKKDPAEF_03230 | 3370941 | AATTAT | AGTTAT | 50.9% | p201RV |
| JKKDPAEF_03230 | hypothetical protein | JKKDPAEF_03230 | 3371376 | GGGGGTAGATGTTTTACCATAACATTATGAAAGATAAGGATT | GGAGGAGTAAAATTTGAGAATGCATTGGGTATGAATAAAGAT | 67.5% | p210B |
| kpsT | Polysialic acid transport ATP-binding protein KpsT | JKKDPAEF_03235 | 3377950 | TAAGTGGGCG | CAATTGAGCG | 94.7% | p261P |
| kpsM | Polysialic acid transport protein KpsM | JKKDPAEF_03236 | 3379769 | ACTCA | ACTCG | 61.1% | p266F |
| epsM | Type II secretion system protein M | JKKDPAEF_03237 | 3379933 | CAGC | CAGT | 93.4% | p261P, p204B |
| JKKDPAEF_03238 | hypothetical protein | JKKDPAEF_03238 | 3381988 | G | A | 1% | p261P |
| flu_2 | Antigen 43 | JKKDPAEF_04005 | 4191843 | AGACTCC | AGGCTCC | 81.2% | p261P |
| cbtA_3 | Cytoskeleton-binding toxin CbtA | JKKDPAEF_04012 | 4196206 | AGAAT | AGAAC | 7.3% | p069PL |
| lpd | Dihydrolipoyl dehydrogenase | JKKDPAEF_04443 | 4655938 | ATGTGGCAGTTATGGGGGGGGGGCCAGGT | ATGTGGCAGTTATGGGGGGGGGCCAGGT | 4.4% | pBB2A |
| cadC | Transcriptional activator CadC | JKKDPAEF_04546 | 4776422 | CATTGTCTACC | CATGGTCTACC | 7.3% | p331F |
| JKKDPAEF_04704 | hypothetical protein | JKKDPAEF_04704 | 4935216 | T | G | 8% | p186F |
| JKKDPAEF_04747 | hypothetical protein | JKKDPAEF_04747 | 4975436 | TGGCCCTGA | TGACCCTGA | 0.9% | p069FHC |
| idnT | Gnt-II system L-idonate transporter | JKKDPAEF_04762 | 4993278 | C | T | 0.9% | p186F |
| yjhG | D-xylonate dehydratase YjhG | JKKDPAEF_04763 | 4996064 | G | A | 1.9% | p186P |
| lgoD_2 | L-galactonate-5-dehydrogenase | JKKDPAEF_04828 | 5065107 | AATCTGATTCGTCACGGTGGCACGGTGGT | AATTTGATTCGTCACGGCGGCACGGTGGT | 1.8% | p572FHCM |
