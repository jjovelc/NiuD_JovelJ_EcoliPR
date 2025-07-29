import pandas as pd
import sys

# Path to your VCF file
vcf_file = sys.argv[1]

# Function to parse the VCF and extract genotypes
def parse_vcf_to_matrix(vcf_file):
    snp_ids = []
    genotypes = []
    samples = None

    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            # Skip headers
            if line.startswith("#CHROM"):
                # Extract sample names
                headers = line.strip().split("\t")
                samples = headers[9:]  # Samples start from the 10th column
            elif not line.startswith("#"):
                fields = line.strip().split("\t")
                snp_id = f"{fields[0]}_{fields[1]}"  # Combine CHROM and POS for SNP ID
                snp_ids.append(snp_id)

                # Extract genotypes
                raw_genotypes = fields[9:]
                numeric_genotypes = [
                    # Map genotypes to numeric labels (e.g., 0/0 -> 0, 0/1 -> 1, 1/1 -> 2, etc.)
                    int(gt.split(":")[0].replace(".", "-1").replace("/", "")) 
                    if gt.split(":")[0] != "./." else -1
                    for gt in raw_genotypes
                ]
                genotypes.append(numeric_genotypes)
    
    # Create a DataFrame
    genotype_matrix = pd.DataFrame(genotypes, index=snp_ids, columns=samples)
    return genotype_matrix

# Generate the matrix
genotype_matrix = parse_vcf_to_matrix(vcf_file)

# Save to a file
genotype_matrix.to_csv("genotype_matrix.csv", sep="\t", na_rep="NA")

print("Genotype matrix created and saved as 'genotype_matrix.csv'")

