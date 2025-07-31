# =============================================================================
# FINAL GENE MAPPING USING SNPEFF GFF FILE
# =============================================================================

setwd('/Users/juanjovel/Library/CloudStorage/OneDrive-UniversityofCalgary/jj/UofC/data_analysis/dongyanNiu/BICOSS')

# Function to extract gene name from snpEff GFF attributes
extract_gene_name_snpeff <- function(attributes) {
  # Split by ";" 
  fields <- strsplit(attributes, ";")[[1]]
  
  # Look for gene= field first (this has the common names like thrA, thrB, etc.)
  gene_field <- fields[grepl("^gene=", fields)]
  
  if (length(gene_field) > 0) {
    gene_name <- gsub("^gene=", "", gene_field[1])
    return(gene_name)
  }
  
  # If no gene= field, try Name= field
  name_field <- fields[grepl("^Name=", fields)]
  
  if (length(name_field) > 0) {
    gene_name <- gsub("^Name=", "", name_field[1])
    return(gene_name)
  }
  
  # If neither, try locus_tag= field
  locus_field <- fields[grepl("^locus_tag=", fields)]
  
  if (length(locus_field) > 0) {
    locus_tag <- gsub("^locus_tag=", "", locus_field[1])
    return(locus_tag)
  }
  
  return(NA)
}

# Function to handle partial coordinates
parse_coordinate <- function(coord_str) {
  if (is.na(coord_str) || coord_str == "" || coord_str == ".") {
    return(NA)
  }
  clean_coord <- gsub("[<>]", "", coord_str)
  result <- suppressWarnings(as.numeric(clean_coord))
  return(result)
}

# Read BICOSS results
cat("=== READING BICOSS RESULTS ===\n")
results <- read.table("BICOSS_all_traits_combined.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Extract SNP positions
results$SNP_Position <- as.numeric(sapply(strsplit(results$SNP_ID, "_"), function(x) x[length(x)]))

cat("Read", nrow(results), "SNPs\n")

# Read the snpEff GFF file
cat("=== READING SNPEFF GFF FILE ===\n")
gff_file <- "NZ_HG941718.1_genes.gff"

if (!file.exists(gff_file)) {
  stop("snpEff GFF file not found: ", gff_file, 
       "\nMake sure you have the snpEff data directory with genes.gff")
}

gff_lines <- readLines(gff_file)
gff_data <- gff_lines[!grepl("^#", gff_lines) & gff_lines != ""]

cat("Read", length(gff_data), "lines from snpEff GFF file\n")

# Initialize result columns
results$Gene_Name <- NA
results$Strand <- NA
results$Gene_Start <- NA
results$Gene_End <- NA
results$Distance_To_Gene <- NA
results$Locus_Tag <- NA

# Process each SNP
cat("=== MAPPING SNPS TO GENES (using snpEff coordinates) ===\n")

for (i in 1:nrow(results)) {
  snp_pos <- results$SNP_Position[i]
  
  if (is.na(snp_pos)) {
    next
  }
  
  # First pass: look for CDS features that CONTAIN the SNP
  # (snpEff GFF uses CDS, not gene features)
  containing_gene <- NULL
  
  for (gff_line in gff_data) {
    fields <- strsplit(gff_line, "\t")[[1]]
    
    if (length(fields) >= 9 && fields[3] == "CDS") {
      start_pos <- parse_coordinate(fields[4])
      end_pos <- parse_coordinate(fields[5])
      
      if (!is.na(start_pos) && !is.na(end_pos)) {
        # Check if SNP is within this CDS
        if (snp_pos >= start_pos && snp_pos <= end_pos) {
          gene_name <- extract_gene_name_snpeff(fields[9])
          
          if (!is.na(gene_name)) {
            # Also extract locus_tag for reference
            locus_match <- regmatches(fields[9], regexpr("locus_tag=([^;]+)", fields[9], perl = TRUE))
            locus_tag <- if (length(locus_match) > 0) gsub("locus_tag=", "", locus_match) else NA
            
            containing_gene <- list(
              name = gene_name,
              strand = fields[7],
              start = start_pos,
              end = end_pos,
              distance = 0,
              locus_tag = locus_tag
            )
            break  # Found containing gene, stop looking
          }
        }
      }
    }
  }
  
  # If we found a containing gene, use it
  if (!is.null(containing_gene)) {
    results$Gene_Name[i] <- containing_gene$name
    results$Strand[i] <- containing_gene$strand
    results$Gene_Start[i] <- containing_gene$start
    results$Gene_End[i] <- containing_gene$end
    results$Distance_To_Gene[i] <- 0
    results$Locus_Tag[i] <- containing_gene$locus_tag
  } else {
    # Second pass: look for closest gene within 2kb
    closest_gene <- NULL
    min_distance <- 2000  # 2kb search window
    
    for (gff_line in gff_data) {
      fields <- strsplit(gff_line, "\t")[[1]]
      
      if (length(fields) >= 9 && fields[3] == "CDS") {
        start_pos <- parse_coordinate(fields[4])
        end_pos <- parse_coordinate(fields[5])
        
        if (!is.na(start_pos) && !is.na(end_pos)) {
          # Calculate distance to gene
          distance <- if (snp_pos < start_pos) {
            start_pos - snp_pos  # upstream
          } else {
            snp_pos - end_pos  # downstream
          }
          
          if (distance <= min_distance) {
            gene_name <- extract_gene_name_snpeff(fields[9])
            
            if (!is.na(gene_name)) {
              locus_match <- regmatches(fields[9], regexpr("locus_tag=([^;]+)", fields[9], perl = TRUE))
              locus_tag <- if (length(locus_match) > 0) gsub("locus_tag=", "", locus_match) else NA
              
              closest_gene <- list(
                name = gene_name,
                strand = fields[7],
                start = start_pos,
                end = end_pos,
                distance = distance,
                locus_tag = locus_tag
              )
              min_distance <- distance
            }
          }
        }
      }
    }
    
    # Use closest gene if found
    if (!is.null(closest_gene)) {
      results$Gene_Name[i] <- closest_gene$name
      results$Strand[i] <- closest_gene$strand
      results$Gene_Start[i] <- closest_gene$start
      results$Gene_End[i] <- closest_gene$end
      results$Distance_To_Gene[i] <- closest_gene$distance
      results$Locus_Tag[i] <- closest_gene$locus_tag
    }
  }
  
  if (i %% 10 == 0) {
    cat("Processed", i, "of", nrow(results), "SNPs\n")
  }
}

# Add copy numbers for duplicated gene names
add_copy_numbers <- function(gene_names) {
  gene_counts <- table(gene_names[!is.na(gene_names)])
  copy_counter <- list()
  result <- character(length(gene_names))
  
  for (i in seq_along(gene_names)) {
    gene <- gene_names[i]
    
    if (is.na(gene)) {
      result[i] <- NA
    } else if (gene_counts[gene] > 1) {
      if (is.null(copy_counter[[gene]])) {
        copy_counter[[gene]] <- 1
      } else {
        copy_counter[[gene]] <- copy_counter[[gene]] + 1
      }
      result[i] <- paste0(gene, "_", copy_counter[[gene]])
    } else {
      result[i] <- gene
    }
  }
  
  return(result)
}

results$Gene_Name_With_Copy <- add_copy_numbers(results$Gene_Name)

# Save results
output_file <- "BICOSS_results_FINAL_gene_names.txt"
write.table(results, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Show results summary
cat("\n=== FINAL RESULTS SUMMARY ===\n")
cat("Results saved to:", output_file, "\n")

mapped_snps <- sum(!is.na(results$Gene_Name))
cat("Total SNPs:", nrow(results), "\n")
cat("Mapped SNPs:", mapped_snps, "\n")
cat("Mapping success rate:", round(100 * mapped_snps / nrow(results), 1), "%\n")

# Check our specific problem SNP
test_snp <- results[results$SNP_Position == 3370198, ]
if (nrow(test_snp) > 0) {
  cat("\nTEST SNP at position 3370198:\n")
  print(test_snp[, c("Trait", "SNP_Position", "Gene_Name", "Gene_Start", "Gene_End", "Distance_To_Gene", "Locus_Tag")])
}

# Show sample results
cat("\nSample results with gene names:\n")
sample_cols <- c("Trait", "SNP_Position", "Gene_Name", "Gene_Name_With_Copy", "Strand", "Distance_To_Gene")
mapped_results <- results[!is.na(results$Gene_Name), ]

if (nrow(mapped_results) > 0) {
  print(head(mapped_results[, sample_cols], 15))
}

# Show unique gene names found
unique_genes <- unique(results$Gene_Name_With_Copy[!is.na(results$Gene_Name_With_Copy)])
cat("\nUnique gene names found:", length(unique_genes), "\n")
cat("Sample gene names (should include thrA, thrB, thrC, etc.):\n")
print(head(sort(unique_genes), 20))

# Create gene summary
if (sum(!is.na(results$Gene_Name)) > 0) {
  gene_summary <- table(results$Gene_Name_With_Copy[!is.na(results$Gene_Name_With_Copy)])
  gene_summary_df <- data.frame(
    Gene = names(gene_summary),
    Count = as.numeric(gene_summary),
    stringsAsFactors = FALSE
  )
  gene_summary_df <- gene_summary_df[order(gene_summary_df$Count, decreasing = TRUE), ]
  
  cat("\nTop genes by association count:\n")
  print(head(gene_summary_df, 10))
  
  write.table(gene_summary_df, "gene_association_summary_FINAL.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Gene summary saved to: gene_association_summary_FINAL.txt\n")
}

cat("\n=== FINAL GENE MAPPING COMPLETE ===\n")
cat("Now you should have proper gene names like thrA, thrB, thrC, yaaA, etc.!\n")
cat("These match exactly what snpEff used for annotation.\n")