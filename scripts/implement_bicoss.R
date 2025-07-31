library(GWAS.BAYES)

# Set working directory
setwd('/Users/juanjovel/Library/CloudStorage/OneDrive-UniversityofCalgary/jj/UofC/data_analysis/dongyanNiu/BICOSS')

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

# Enhanced data validation function
validate_bicoss_inputs <- function(Y, SNPs, kinship) {
  cat("=== Input Validation ===\n")
  
  # Check Y
  cat("Y - Length:", length(Y), "Type:", class(Y), "\n")
  cat("Y - Range:", range(Y, na.rm = TRUE), "\n")
  cat("Y - NA count:", sum(is.na(Y)), "\n")
  cat("Y - Variance:", var(Y, na.rm = TRUE), "\n")
  
  # Check SNPs
  cat("SNPs - Dimensions:", dim(SNPs), "Type:", class(SNPs), "\n")
  cat("SNPs - Range:", range(SNPs, na.rm = TRUE), "\n")
  cat("SNPs - NA count:", sum(is.na(SNPs)), "\n")
  
  # Check kinship
  cat("Kinship - Dimensions:", dim(kinship), "Type:", class(kinship), "\n")
  cat("Kinship - Range:", range(kinship, na.rm = TRUE), "\n")
  cat("Kinship - NA count:", sum(is.na(kinship)), "\n")
  cat("Kinship - Is symmetric:", isSymmetric(kinship), "\n")
  
  # Check dimension compatibility
  cat("Dimension checks:\n")
  cat("  Y length == SNPs rows:", length(Y) == nrow(SNPs), "\n")
  cat("  Y length == kinship rows:", length(Y) == nrow(kinship), "\n")
  cat("  Y length == kinship cols:", length(Y) == ncol(kinship), "\n")
  
  return(TRUE)
}

# Main loop with enhanced error handling
multivariate_results <- list()

for (trait in colnames(traits_data)) {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("Processing trait:", trait, "\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  tryCatch({
    # Prepare phenotype data
    Y <- as.numeric(as.character(traits_data[[trait]]))
    valid_samples <- !is.na(Y)
    Y <- Y[valid_samples]
    
    # Check if we have enough valid samples
    if (sum(valid_samples) < 10) {
      cat("Warning: Too few valid samples for trait", trait, ". Skipping.\n")
      next
    }
    
    # Center Y but don't scale to unit variance yet
    Y <- as.numeric(scale(Y, center = TRUE, scale = FALSE))
    
    # Prepare genotype data
    geno_matrix_valid <- geno_matrix[valid_samples, , drop = FALSE]
    
    # More stringent SNP filtering
    cat("Filtering SNPs...\n")
    snp_var <- apply(geno_matrix_valid, 2, function(x) {
      if (all(is.na(x))) return(0)
      return(var(x, na.rm = TRUE))
    })
    
    # Remove SNPs with zero variance or too many missing values
    missing_prop <- apply(geno_matrix_valid, 2, function(x) sum(is.na(x)) / length(x))
    valid_snps <- !is.na(snp_var) & snp_var > 1e-6 & missing_prop < 0.1
    
    cat("SNPs before filtering:", ncol(geno_matrix_valid), "\n")
    cat("SNPs after variance filtering:", sum(valid_snps), "\n")
    
    if (sum(valid_snps) < 10) {
      cat("Warning: Too few valid SNPs for trait", trait, ". Skipping.\n")
      next
    }
    
    # Keep track of original SNP names
    snp_names <- colnames(geno_matrix_valid)[valid_snps]
    
    # Subset and properly scale genotype matrix
    geno_matrix_valid <- geno_matrix_valid[, valid_snps, drop = FALSE]
    
    # Handle any remaining NAs
    if (any(is.na(geno_matrix_valid))) {
      cat("Warning: NAs found in genotype matrix. Imputing with column means.\n")
      for (j in 1:ncol(geno_matrix_valid)) {
        col_mean <- mean(geno_matrix_valid[, j], na.rm = TRUE)
        geno_matrix_valid[is.na(geno_matrix_valid[, j]), j] <- col_mean
      }
    }
    
    # Scale genotype matrix
    geno_matrix_valid <- scale(geno_matrix_valid, center = TRUE, scale = TRUE)
    
    # Prepare kinship matrix
    kinship_matrix_valid <- kinship_matrix[valid_samples, valid_samples, drop = FALSE]
    
    # Ensure kinship matrix is positive definite
    eigen_vals <- eigen(kinship_matrix_valid, symmetric = TRUE, only.values = TRUE)$values
    min_eigen <- min(eigen_vals)
    if (min_eigen <= 0) {
      cat("Warning: Kinship matrix not positive definite. Adding regularization.\n")
      diag(kinship_matrix_valid) <- diag(kinship_matrix_valid) + abs(min_eigen) + 1e-6
    }
    
    # Validate inputs before BICOSS call
    validate_bicoss_inputs(Y, geno_matrix_valid, kinship_matrix_valid)
    
    # Convert to matrices and ensure proper types
    Y_vec <- as.numeric(Y)
    SNPs_mat <- as.matrix(geno_matrix_valid)
    kinship_mat <- as.matrix(kinship_matrix_valid)
    
    # Ensure no infinite values
    if (any(!is.finite(Y_vec))) {
      cat("Error: Non-finite values in Y. Skipping trait.\n")
      next
    }
    if (any(!is.finite(SNPs_mat))) {
      cat("Error: Non-finite values in SNPs matrix. Skipping trait.\n")
      next
    }
    if (any(!is.finite(kinship_mat))) {
      cat("Error: Non-finite values in kinship matrix. Skipping trait.\n")
      next
    }
    
    cat("Calling BICOSS function...\n")
    set.seed(123)
    
    # Try with more conservative parameters first
    BICOSS_result <- BICOSS(
      Y = Y_vec,
      SNPs = SNPs_mat,
      kinship = kinship_mat,
      FDR_Nominal = 0.1,  # Less stringent initially
      P3D = TRUE,         # Try with P3D = TRUE first
      maxiterations = 200,  # Fewer iterations initially
      runs_til_stop = 20    # Stop earlier
    )
    
    cat("BICOSS completed successfully!\n")
    cat("Best model length:", length(BICOSS_result$best_model), "\n")
    
    # Process results
    if (!is.null(BICOSS_result$best_model) && 
        length(BICOSS_result$best_model) > 0 && 
        !all(is.na(BICOSS_result$best_model))) {
      
      cat("Found significant SNPs:", paste(BICOSS_result$best_model, collapse=", "), "\n")
      
      # Calculate effects for significant SNPs
      if (length(BICOSS_result$best_model) == 1) {
        snp_data <- SNPs_mat[, BICOSS_result$best_model, drop=FALSE]
        cor_test <- cor.test(as.numeric(snp_data), Y_vec, method="pearson")
        
        significant_snps <- data.frame(
          Trait = trait,
          SNP_Index = BICOSS_result$best_model,
          SNP_ID = snp_names[BICOSS_result$best_model],
          Effect = cor_test$estimate,
          P_Value = cor_test$p.value,
          stringsAsFactors = FALSE
        )
      } else {
        effects <- numeric(length(BICOSS_result$best_model))
        pvalues <- numeric(length(BICOSS_result$best_model))
        
        for(i in seq_along(BICOSS_result$best_model)) {
          snp_idx <- BICOSS_result$best_model[i]
          snp_data <- SNPs_mat[, snp_idx, drop=FALSE]
          cor_test <- cor.test(as.numeric(snp_data), Y_vec, method="pearson")
          effects[i] <- cor_test$estimate
          pvalues[i] <- cor_test$p.value
        }
        
        significant_snps <- data.frame(
          Trait = trait,
          SNP_Index = BICOSS_result$best_model,
          SNP_ID = snp_names[BICOSS_result$best_model],
          Effect = effects,
          P_Value = pvalues,
          stringsAsFactors = FALSE
        )
      }
      
      # Save results
      outfile <- paste0("BICOSS_", trait, "_significant_SNPs.txt")
      write.table(significant_snps, file = outfile, 
                  sep = "\t", quote = FALSE, row.names = FALSE)
      
      cat("Results saved to:", outfile, "\n")
      print(significant_snps)
      
      # Store in multivariate results
      multivariate_results[[trait]] <- significant_snps
      
    } else {
      cat("No significant SNPs found for trait:", trait, "\n")
      
      # Create empty results file
      outfile <- paste0("BICOSS_", trait, "_significant_SNPs.txt")
      empty_df <- data.frame(
        Trait = character(0),
        SNP_Index = numeric(0),
        SNP_ID = character(0),
        Effect = numeric(0),
        P_Value = numeric(0)
      )
      write.table(empty_df, file = outfile, 
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }
    
  }, error = function(e) {
    cat("\n!!! ERROR processing trait:", trait, "!!!\n")
    cat("Error message:", conditionMessage(e), "\n")
    cat("Error class:", class(e), "\n")
    
    # Print more detailed error information
    if (exists("traceback")) {
      cat("Traceback:\n")
      print(traceback())
    }
    
    # Create empty results file for failed trait
    outfile <- paste0("BICOSS_", trait, "_significant_SNPs.txt")
    empty_df <- data.frame(
      Trait = character(0),
      SNP_Index = numeric(0),
      SNP_ID = character(0),
      Effect = numeric(0),
      P_Value = numeric(0)
    )
    write.table(empty_df, file = outfile, 
                sep = "\t", quote = FALSE, row.names = FALSE)
  })
}

# Save combined results
if (length(multivariate_results) > 0) {
  cat("\n", paste(rep("=", 50), collapse=""), "\n")
  cat("SAVING COMBINED RESULTS\n")
  cat(paste(rep("=", 50), collapse=""), "\n")
  
  # Combine all results
  combined_results <- do.call(rbind, multivariate_results)
  write.table(combined_results, file = "BICOSS_all_traits_significant_SNPs.txt", 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("Combined results saved to: BICOSS_all_traits_significant_SNPs.txt\n")
  cat("Total traits with results:", length(multivariate_results), "\n")
  cat("Total significant SNPs:", nrow(combined_results), "\n")
} else {
  cat("\nNo successful results to combine.\n")
}

##### ANALYSIS COMPLETE
cat("\n")
cat(paste(rep("=", 50), collapse=""), "\n")
cat("BICOSS ANALYSIS COMPLETED\n")
cat(paste(rep("=", 50), collapse=""), "\n")
cat("Results summary:\n")
cat("- Processed", length(colnames(traits_data)), "traits\n")
cat("- Successful results:", length(multivariate_results), "traits\n")
cat("- Individual files: BICOSS_[trait]_significant_SNPs.txt\n")
if (length(multivariate_results) > 0) {
  cat("- Combined file: BICOSS_all_traits_significant_SNPs.txt\n")
}
cat("- Kinship matrix: kinship_matrix.txt\n")
cat(paste(rep("=", 50), collapse=""), "\n")
