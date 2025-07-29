library(data.table)

# Load the genotype data
geno_data <- fread("genotypes.tsv", header = FALSE)
colnames(geno_data) <- c("SNP", "Sample", "Genotype")

# Convert genotypes to numeric values
geno_data$Genotype <- as.numeric(gsub("0/0", "0", 
                               gsub("0/1", "1", 
                               gsub("1/1", "2", 
                               gsub("\\./\\.", "NA", geno_data$Genotype)))))

