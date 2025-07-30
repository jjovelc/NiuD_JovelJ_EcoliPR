library(GWAS.BAYES)

# Set working directory
setwd('/Users/juanjovel/OneDrive/jj/UofC/data_analysis/dongyanNiu/BICOSS')

# Define the input file
infile <- "genotype_matrix.tsv"
traits_file <- "phage_traits.tsv"

# Read the data
indata <- read.table(infile, sep = "\t", row.names = 1, header = TRUE)
traits_data <- read.table(traits_file, sep = "\t", row.names = 1, header = TRUE)

# Print initial dimensions to verify
cat("Initial dimensions:\n")
cat("indata (SNPs x Samples):", dim(indata), "\n")
cat("traits_data:", dim(traits_data), "\n")

# Important: Transpose indata to get samples as rows
geno_matrix <- t(indata)
cat("geno_matrix (Samples x SNPs):", dim(geno_matrix), "\n")

# Scale the genotype matrix
geno_matrix <- scale(geno_matrix, center = TRUE, scale = FALSE)

# Calculate kinship matrix correctly (samples x samples)
kinship_matrix <- tcrossprod(geno_matrix) / ncol(geno_matrix)
cat("kinship_matrix (Samples x Samples):", dim(kinship_matrix), "\n")

# Verify dimensions
stopifnot(
  nrow(kinship_matrix) == nrow(geno_matrix),
  ncol(kinship_matrix) == nrow(geno_matrix)
)

# Write the correctly sized kinship matrix
write.table(kinship_matrix, file = "kinship_matrix.txt", 
            sep = "\t", quote = FALSE, 
            row.names = TRUE, col.names = TRUE)


# Now the loop to process each phage
# Main loop with proper SNP filtering
for (trait in colnames(traits_data)) {
  cat("\n\nProcessing trait:", trait, "\n")
  
  # Prepare phenotype data
  Y <- as.numeric(as.character(traits_data[[trait]]))
  valid_samples <- !is.na(Y)
  Y <- Y[valid_samples]
  Y <- scale(Y, center = TRUE, scale = FALSE)
  
  # Print Y diagnostics
  cat("Y summary:", "\n")
  print(summary(Y))
  
  # Prepare genotype data
  geno_matrix_valid <- geno_matrix[valid_samples, ]
  
  # Calculate SNP variance
  snp_var <- apply(geno_matrix_valid, 2, var)
  valid_snps <- !is.na(snp_var) & snp_var > 0
  cat("Number of valid SNPs:", sum(valid_snps), "\n")
  
  # Keep track of original SNP names/positions before subsetting
  snp_names <- colnames(geno_matrix_valid)[valid_snps]
  
  # Now subset and scale the genotype matrix
  geno_matrix_valid <- geno_matrix_valid[, valid_snps]
  geno_matrix_valid <- scale(geno_matrix_valid, center = TRUE, scale = TRUE)
  
  # Prepare kinship matrix
  kinship_matrix_valid <- kinship_matrix[valid_samples, valid_samples]
  
  tryCatch({
    set.seed(123)
    
    BICOSS_Exact <- BICOSS(
      Y = as.numeric(Y),
      SNPs = as.matrix(geno_matrix_valid),
      kinship = as.matrix(kinship_matrix_valid),
      FDR_Nominal = 0.05,
      P3D = FALSE,
      maxiterations = 400,
      runs_til_stop = 40
    )
    
    cat("BICOSS run completed\n")
    cat("best_model length:", length(BICOSS_Exact$best_model), "\n")
    cat("best_model content:", paste(BICOSS_Exact$best_model, collapse=", "), "\n")
    
    # Check if best_model contains valid SNP indices
    if (!is.null(BICOSS_Exact$best_model) && 
        length(BICOSS_Exact$best_model) > 0 && 
        !is.character(BICOSS_Exact$best_model)) {
      
      # Handle single SNP case
      if (length(BICOSS_Exact$best_model) == 1) {
        snp_data <- geno_matrix_valid[, BICOSS_Exact$best_model, drop=FALSE]
        cor_test <- cor.test(snp_data, Y, method="pearson")
        
        significant_snps <- data.frame(
          Trait = trait,
          SNP_Index = BICOSS_Exact$best_model,
          SNP_ID = snp_names[BICOSS_Exact$best_model],
          Effect = cor_test$estimate,
          P_Value = cor_test$p.value
        )
      } else {
        # Handle multiple SNPs case
        effects <- numeric(length(BICOSS_Exact$best_model))
        pvalues <- numeric(length(BICOSS_Exact$best_model))
        
        for(i in seq_along(BICOSS_Exact$best_model)) {
          snp_data <- geno_matrix_valid[, BICOSS_Exact$best_model[i], drop=FALSE]
          cor_test <- cor.test(snp_data, Y, method="pearson")
          effects[i] <- cor_test$estimate
          pvalues[i] <- cor_test$p.value
        }
        
        significant_snps <- data.frame(
          Trait = trait,
          SNP_Index = BICOSS_Exact$best_model,
          SNP_ID = snp_names[BICOSS_Exact$best_model],
          Effect = effects,
          P_Value = pvalues
        )
      }
      
      # Save results
      outfile <- paste0("BICOSS_", trait, "_significant_SNPs.txt")
      write.table(significant_snps, file = outfile, 
                  sep = "\t", quote = FALSE, row.names = FALSE)
      
      cat("\nResults for trait:", trait, "\n")
      cat("Number of significant SNPs found:", nrow(significant_snps), "\n")
      print(significant_snps)
    } else {
      cat("\nNo significant SNPs found for trait:", trait, "\n")
      outfile <- paste0("BICOSS_", trait, "_significant_SNPs.txt")
      write.table(data.frame(), file = outfile, 
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }, error = function(e) {
    cat("\nError processing trait:", trait, "\n")
    cat("Error message:", conditionMessage(e), "\n")
    cat("Error trace:\n")
    print(e)
    
    outfile <- paste0("BICOSS_", trait, "_significant_SNPs.txt")
    write.table(data.frame(), file = outfile, 
                sep = "\t", quote = FALSE, row.names = FALSE)
  })
}


##### MULTIVARIATE ANALYSIS
# Combine multiple traits into a matrix and pass them as Y.

# Create a multivariate trait matrix
Y_multivariate <- as.matrix(as.numeric(traits_data))

# Run BICOSS analysis
BICOSS_Exact_multi <- BICOSS(
  Y = Y_multivariate,
  SNPs = geno_matrix,
  kinship = kinship_matrix,
  FDR_Nominal = 0.05,
  P3D = FALSE,
  maxiterations = 400,
  runs_til_stop = 40
)

# Save results
write.table(BICOSS_Exact_multi$results, file = "BICOSS_multivariate_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)


