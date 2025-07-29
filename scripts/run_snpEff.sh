#!/usr/bin/bash

java -jar ../snpEff.jar build -gff3 -v NZ_HG941718.1

java -jar ../snpEff.jar -v NZ_HG941718.1 ../../E.coli_NZ_HG941718.1.variants_filtered50-30-0.05.vcf > E.coli_NZ_HG941718.1.variants_filtered50-30-0.05_annotated.vcf
