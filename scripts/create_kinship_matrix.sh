# 1. Exclude synonymous mutations

INFILE=$1  # E.coli_NZ_HG941718.1.variants_filtered50-30-0.05_annotated.vcf
OUTFILE=$2 # E.coli_NZ_HG941718.1_nonSynonymousVariants.vcf 

bcftools view -i 'INFO/ANN[0] !~ "synonymous_variant"' "$INFILE" > "$OUTFILE"

# 2. Extract genotypes

INFILE=$1  # E.coli_NZ_HG941718.1_nonSynonymousVariants.vcf
OUTFILE=$2 # E.coli_NZ_HG941718.1_nonSynonymousVariants_genotypes.txt 

bcftools query -f '[%CHROM:%POS\t%SAMPLE\t%GT\n]' "$INFILE" > "$OUTFILE"


