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
        
        return(list(
          position = pos,
          ref_allele = ref_allele,
          alt_allele = alt_allele,
          allele_freq = allele_freq,
          variant_type = variant_type,
          product_name = product_name,
          protein_id = protein_id
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
  Protein_Product = character(0),
  Protein_ID = character(0),
  Reference_Genome_Mutation_Site = character(0),
  Reference_Allele = character(0),
  SNP_Allele = character(0),
  Percentage_of_isolates_carrying_this_SNPs = character(0),
  Traits_Associated = character(0),
  Effect_Size = character(0),
  P_Value = character(0),
  stringsAsFactors = FALSE
)

# Process each significant result
for (i in 1:nrow(significant_results)) {
  
  position <- significant_results$SNP_Position[i]
  gene_name <- significant_results$Gene_Name[i]
  locus_tag <- significant_results$Locus_Tag[i]
  trait <- significant_results$Trait[i]
  effect <- significant_results$Effect[i]
  pvalue <- significant_results$P_Value[i]
  
  cat("Processing position", position, "for gene", gene_name, "\n")
  
  # Get VCF information
  vcf_info <- extract_vcf_info(position, vcf_file)
  
  # Get product name from GFF
  product_name <- get_product_from_gff(locus_tag, gff_file)
  
  if (!is.null(vcf_info)) {
    
    # Use gene name as product if no specific product found
    if (is.na(product_name) || product_name == "") {
      product_name <- gene_name
    }
    
    # Calculate percentage of isolates carrying the SNP
    percentage <- if (!is.na(vcf_info$allele_freq)) {
      paste0(round(vcf_info$allele_freq * 100, 1), "%")
    } else {
      "Unknown"
    }
    
    # Handle multiple alternative alleles
    alt_alleles <- strsplit(vcf_info$alt_allele, ",")[[1]]
    main_alt <- alt_alleles[1]  # Take first alternative allele
    
    # Add row to summary table
    summary_table <- rbind(summary_table, data.frame(
      Protein_Product = product_name,
      Protein_ID = ifelse(is.na(locus_tag), gene_name, locus_tag),
      Reference_Genome_Mutation_Site = vcf_info$position,
      Reference_Allele = vcf_info$ref_allele,
      SNP_Allele = main_alt,
      Percentage_of_isolates_carrying_this_SNPs = percentage,
      Traits_Associated = trait,
      Effect_Size = round(effect, 3),
      P_Value = formatC(pvalue, format = "e", digits = 2),
      stringsAsFactors = FALSE
    ))
  }
}

# Group by position to combine traits for the same SNP
cat("\n=== COMBINING TRAITS FOR SAME SNPs ===\n")

final_table <- summary_table %>%
  group_by(Reference_Genome_Mutation_Site, Reference_Allele, SNP_Allele) %>%
  summarise(
    Protein_Product = first(Protein_Product),
    Protein_ID = first(Protein_ID),
    Percentage_of_isolates_carrying_this_SNPs = first(Percentage_of_isolates_carrying_this_SNPs),
    Traits_Associated = paste(unique(Traits_Associated), collapse = ", "),
    Effect_Size = paste(unique(Effect_Size), collapse = ", "),
    P_Value = paste(unique(P_Value), collapse = ", "),
    .groups = 'drop'
  ) %>%
  select(
    Protein_Product,
    Protein_ID, 
    Reference_Genome_Mutation_Site,
    Reference_Allele,
    SNP_Allele,
    Percentage_of_isolates_carrying_this_SNPs,
    Traits_Associated,
    Effect_Size,
    P_Value
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
  "Protein Product",
  "Protein ID", 
  "Reference Genome Mutation Site",
  "Reference Allele",
  "SNP Allele",
  "Percentage of isolates carrying this SNPs",
  "Traits Associated",
  "Effect Size",
  "P-Value"
)

# Save formatted version
write.table(formatted_table, "Variant_Summary_Table_Formatted.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Also save as CSV for Excel
write.csv(formatted_table, "Variant_Summary_Table.csv", row.names = FALSE)

cat("Formatted table saved to: Variant_Summary_Table_Formatted.txt\n")
cat("CSV version saved to: Variant_Summary_Table.csv\n")

# Show statistics
cat("\n=== TABLE STATISTICS ===\n")
cat("Total variants included:", nrow(final_table), "\n")
cat("Unique genes affected:", length(unique(final_table$Protein_Product)), "\n")
cat("Unique traits involved:", length(unique(unlist(strsplit(final_table$Traits_Associated, ", ")))), "\n")

cat("\n=== VARIANT SUMMARY TABLE GENERATION COMPLETE ===\n")
