eval "$(conda shell.bash hook)"
conda activate call_variants_nano

for FILE in E.coli_NZ_HG941718.1.variants.vcf; do biopet-vcffilter -I ${FILE} -o ${FILE/.vcf/_filtered50-30-0.05.vcf} --minSampleDepth 50 --minQualScore 30 --minAlternateDepth 3; done
