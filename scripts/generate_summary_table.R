# =============================================================================
# GENERATE VARIANT SUMMARY TABLE FROM VCF AND GENE RESULTS
# =============================================================================

# Load required libraries
library(dplyr)

# Set working directory
setwd('/Users/juanjovel/Library/CloudStorage/OneDrive-UniversityofCalgary/jj/UofC/data_analysis/dongyanNiu/BICOSS')

# Read the final BICOSS results with gene names
cat("=== READING BICOSS RESULTS WITH GENE NAMES ===\n")
results <- read.table("BICOSS_results_FINAL_gene_names.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Filter for significant associations with mapped genes
significant_results <- results %>%
  filter(!is.na(Gene_Name) & P_Value < 0.05) %>%
  arrange(P_Value)

cat("Found", nrow(significant_results), "significant associations with gene names\n")

# Function to extract variant information from VCF
extract_vcf_info <- function(position, vcf_file) {
  
  # Search for the position in VCF
  cmd <- paste0("grep '^NZ_HG941718.1\t", position, "\t' ", vcf_file)
  
  tryCatch({
    vcf_line <- system(cmd, intern = TRUE)
    
    if (length(vcf_line) > 0) {
      # Split VCF line into fields
      fields <- strsplit(vcf_line, "\t")[[1]]
      
      if (length(fields) >= 8) {
        chrom <- fields[1]
        pos <- fields[2]
        ref_allele <- fields[4]
        alt_allele <- fields[5]
        info_field <- fields[8]
        
        # Extract allele frequency (AF) from INFO field
        af_match <- regmatches(info_field, regexpr("AF=([^;]+)", info_field, perl = TRUE))
        allele_freq <- if (length(af_match) > 0) {
          af_value <- gsub("AF=", "", af_match)
          # Handle multiple allele frequencies (comma-separated)
          af_values <- strsplit(af_value, ",")[[1]]
          as.numeric(af_values[1])  # Take first allele frequency
        } else {
          NA
        }
        
        # Extract variant type from TYPE field
        type_match <- regmatches(info_field, regexpr("TYPE=([^;]+)", info_field, perl = TRUE))
        variant_type <- if (length(type_match) > 0) {
          gsub("TYPE=", "", type_match)
        } else {
          "unknown"
        }
        
        # Extract annotation from ANN field
        product_name <- NA
        protein_id <- NA
        
        if (grepl("ANN=", info_field)) {
          ann_match <- regmatches(info_field, regexpr("ANN=([^;]+)", info_field, perl = TRUE))
          ann_content <- gsub("ANN=", "", ann_match)
          
          # Split by comma (multiple annotations) and take first one
          first_annotation <- strsplit(ann_content, ",")[[1]][1]
          
          # Split by pipe to get fields
          ann_fields <- strsplit(first_annotation, "\\|")[[1]]
          
          if (length(ann_fields) >= 4) {
            # Try to get product name from various fields
            if (length(ann_fields) >= 11 && ann_fields[11] != "") {
              product_name <- ann_fields[11]  # HGVS.c field sometimes has product info
            }
            
            # Extract locus tag as protein ID
            if (ann_fields[4] != "" && ann_fields[4] != "null") {
              protein_id <- ann_fields[4]
            }
          }
        }
        
        # Calculate percentage of isolates carrying the SNP from genotype data
        # Get the genotype data (format field and sample data)
        if (length(fields) >= 9) {
          format_field <- fields[9]
          sample_data <- fields[10:length(fields)]
          
          # Parse genotype data to count samples carrying the variant
          samples_with_variant <- 0
          total_samples <- 0
          
          for (sample in sample_data) {
            if (sample != "." && sample != "./.") {
              total_samples <- total_samples + 1
              
              # Extract genotype (first field before colon)
              genotype <- strsplit(sample, ":")[[1]][1]
              
              # Check if sample carries the variant (genotype 1 or 2)
              if (grepl("1", genotype) || grepl("2", genotype)) {
                samples_with_variant <- samples_with_variant + 1
              }
            }
          }
          
          # Calculate percentage
          if (total_samples > 0) {
            percentage_carrying <- (samples_with_variant / total_samples) * 100
          } else {
            percentage_carrying <- NA
          }
        } else {
          percentage_carrying <- NA
        }
        
        return(list(
          position = pos,
          ref_allele = ref_allele,
          alt_allele = alt_allele,
          allele_freq = allele_freq,
          variant_type = variant_type,
          product_name = product_name,
          protein_id = protein_id,
          percentage_carrying = percentage_carrying
        ))
      }
    }
  }, error = function(e) {
    return(NULL)
  })
  
  return(NULL)
}

# Function to get product name from GFF using locus tag
get_product_from_gff <- function(locus_tag, gff_file) {
  
  if (is.na(locus_tag) || locus_tag == "") {
    return(NA)
  }
  
  tryCatch({
    # Search for the locus tag in GFF
    cmd <- paste0("grep '", locus_tag, "' ", gff_file)
    gff_lines <- system(cmd, intern = TRUE)
    
    if (length(gff_lines) > 0) {
      # Take first matching line
      gff_line <- gff_lines[1]
      fields <- strsplit(gff_line, "\t")[[1]]
      
      if (length(fields) >= 9) {
        attributes <- fields[9]
        
        # Extract product information
        product_match <- regmatches(attributes, regexpr("product=([^;]+)", attributes, perl = TRUE))
        if (length(product_match) > 0) {
          product <- gsub("product=", "", product_match)
          return(product)
        }
      }
    }
  }, error = function(e) {
    return(NA)
  })
  
  return(NA)
}

# Create the summary table
cat("\n=== CREATING VARIANT SUMMARY TABLE ===\n")

vcf_file <- "E.coli_NZ_HG941718.1.variants_filtered50-30-0.05_annotated.vcf"
gff_file <- "NZ_HG941718.1_genes.gff"

# Initialize the summary table
summary_table <- data.frame(
  Gene_name = character(0),
  Protein_Product = character(0),
  Protein_ID = character(0),
  Reference_Genome_Mutation_Site = character(0),
  Reference_Allele = character(0),
  SNP_Allele = character(0),
  Percentage_of_isolates_carrying_this_SNPs = character(0),
  Traits_Associated = character(0),
  stringsAsFactors = FALSE
)

# Process each significant result
for (i in 1:nrow(significant_results)) {
  
  position <- significant_results$SNP_Position[i]
  gene_name <- significant_results$Gene_Name[i]
  locus_tag <- significant_results$Locus_Tag[i]
  trait <- significant_results$Trait[i]
  
  cat("Processing position", position, "for gene", gene_name, "\n")
  
  # Get VCF information
  vcf_info <- extract_vcf_info(position, vcf_file)
  
  # Get product name from GFF
  product_name <- get_product_from_gff(locus_tag, gff_file)
  
  if (!is.null(vcf_info)) {
    
    # Include all variants regardless of frequency
    # (Previously filtered out high-frequency variants > 50%)
    cat("  Including variant with frequency:", round(vcf_info$percentage_carrying, 1), "%\n")
    
    # Use gene name as product if no specific product found
    if (is.na(product_name) || product_name == "") {
      product_name <- gene_name
    }
    
    # Calculate percentage of isolates carrying the SNP
    percentage <- if (!is.na(vcf_info$percentage_carrying)) {
      paste0(round(vcf_info$percentage_carrying, 1), "%")
    } else {
      "Unknown"
    }
    
    # Handle multiple alternative alleles
    alt_alleles <- strsplit(vcf_info$alt_allele, ",")[[1]]
    main_alt <- alt_alleles[1]  # Take first alternative allele
    
    # Add row to summary table
    summary_table <- rbind(summary_table, data.frame(
      Gene_name = gene_name,
      Protein_Product = product_name,
      Protein_ID = ifelse(is.na(locus_tag), gene_name, locus_tag),
      Reference_Genome_Mutation_Site = vcf_info$position,
      Reference_Allele = vcf_info$ref_allele,
      SNP_Allele = main_alt,
      Percentage_of_isolates_carrying_this_SNPs = percentage,
      Traits_Associated = trait,
      stringsAsFactors = FALSE
    ))
  }
}

# Group by position to combine traits for the same SNP
cat("\n=== COMBINING TRAITS FOR SAME SNPs ===\n")

final_table <- summary_table %>%
  group_by(Reference_Genome_Mutation_Site, Reference_Allele, SNP_Allele) %>%
  summarise(
    Gene_name = first(Gene_name),
    Protein_Product = first(Protein_Product),
    Protein_ID = first(Protein_ID),
    Percentage_of_isolates_carrying_this_SNPs = first(Percentage_of_isolates_carrying_this_SNPs),
    Traits_Associated = paste(unique(Traits_Associated), collapse = ", "),
    .groups = 'drop'
  ) %>%
  select(
    Gene_name,
    Protein_Product,
    Protein_ID, 
    Reference_Genome_Mutation_Site,
    Reference_Allele,
    SNP_Allele,
    Percentage_of_isolates_carrying_this_SNPs,
    Traits_Associated
  ) %>%
  arrange(as.numeric(Reference_Genome_Mutation_Site))

# Save the summary table
output_file <- "Variant_Summary_Table.txt"
write.table(final_table, output_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== SUMMARY TABLE CREATED ===\n")
cat("Summary table saved to:", output_file, "\n")
cat("Number of unique variants:", nrow(final_table), "\n")

# Display the table
cat("\nVariant Summary Table:\n")
print(final_table)

# Create a formatted version for publication
cat("\n=== CREATING FORMATTED VERSION ===\n")

# Create a nicely formatted version
formatted_table <- final_table
colnames(formatted_table) <- c(
  "Gene Name",
  "Protein Product",
  "Protein ID", 
  "Reference Genome Mutation Site",
  "Reference Allele",
  "SNP Allele",
  "Percentage of isolates carrying this SNPs",
  "Traits Associated"
)

# Save formatted version
write.table(formatted_table, "Variant_Summary_Table_Formatted.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Also save as CSV for Excel
write.csv(formatted_table, "Variant_Summary_Table.csv", row.names = FALSE)

cat("Formatted table saved to: Variant_Summary_Table_Formatted.txt\n")
cat("CSV version saved to: Variant_Summary_Table.csv\n")

# Create markdown formatted table
cat("\n=== CREATING MARKDOWN FORMATTED TABLE ===\n")

# Function to create markdown table
create_markdown_table <- function(df) {
  # Create header
  header <- paste("|", paste(colnames(df), collapse = " | "), "|")
  
  # Create separator
  separator <- paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|")
  
  # Create rows
  rows <- apply(df, 1, function(row) {
    paste("|", paste(row, collapse = " | "), "|")
  })
  
  # Combine all parts
  markdown_table <- c(header, separator, rows)
  
  return(markdown_table)
}

# Create markdown table
markdown_table <- create_markdown_table(formatted_table)

# Save markdown table to file
writeLines(markdown_table, "Variant_Summary_Table_Markdown.md")

cat("Markdown table saved to: Variant_Summary_Table_Markdown.md\n")

# Print markdown table to console
cat("\n=== MARKDOWN FORMATTED TABLE ===\n")
cat(paste(markdown_table, collapse = "\n"))
cat("\n")

# Show statistics
cat("\n=== TABLE STATISTICS ===\n")
cat("Total variants included:", nrow(final_table), "\n")
cat("Unique genes affected:", length(unique(final_table$Gene_name)), "\n")
cat("Unique traits involved:", length(unique(unlist(strsplit(final_table$Traits_Associated, ", ")))), "\n")

# Show frequency distribution
cat("\n=== FREQUENCY DISTRIBUTION ===\n")
if (nrow(final_table) > 0) {
  # Extract percentages and convert to numeric
  percentages <- as.numeric(gsub("%", "", final_table$Percentage_of_isolates_carrying_this_SNPs))
  percentages <- percentages[!is.na(percentages)]
  
  if (length(percentages) > 0) {
    cat("Mean frequency:", round(mean(percentages), 1), "%\n")
    cat("Median frequency:", round(median(percentages), 1), "%\n")
    cat("Min frequency:", round(min(percentages), 1), "%\n")
    cat("Max frequency:", round(max(percentages), 1), "%\n")
    cat("Variants with frequency < 10%:", sum(percentages < 10), "\n")
    cat("Variants with frequency 10-25%:", sum(percentages >= 10 & percentages <= 25), "\n")
    cat("Variants with frequency > 25%:", sum(percentages > 25), "\n")
  }
}

cat("\n=== VARIANT SUMMARY TABLE GENERATION COMPLETE ===\n")