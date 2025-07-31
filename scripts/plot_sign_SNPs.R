# =============================================================================
# CREATE SNP-TRAIT ASSOCIATIONS DOT PLOT
# =============================================================================

# Load required libraries
library(ggplot2)
library(dplyr)

# Set working directory
setwd('/Users/juanjovel/Library/CloudStorage/OneDrive-UniversityofCalgary/jj/UofC/data_analysis/dongyanNiu/BICOSS')

# Read the final results
cat("=== READING FINAL BICOSS RESULTS ===\n")
results <- read.table("BICOSS_results_FINAL_gene_names.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Filter for significant associations and mapped genes
cat("Filtering data for plotting...\n")
plot_data <- results %>%
  filter(!is.na(Gene_Name) & !is.na(P_Value) & P_Value < 0.05) %>%  # Significant associations only
  mutate(
    log10_pvalue = -log10(P_Value),  # Convert p-values to -log10 scale
    abs_effect = abs(Effect),        # Absolute effect size for dot size
    effect_direction = ifelse(Effect > 0, "Positive", "Negative")  # Effect direction for color
  )

cat("Found", nrow(plot_data), "significant associations for plotting\n")

# Check if we have data to plot
if (nrow(plot_data) == 0) {
  stop("No significant associations found to plot. Try relaxing the p-value threshold.")
}

# Show summary of data
cat("\nData summary:\n")
cat("P-value range:", min(plot_data$P_Value), "to", max(plot_data$P_Value), "\n")
cat("Effect size range:", min(plot_data$Effect), "to", max(plot_data$Effect), "\n")
cat("Unique genes:", length(unique(plot_data$Gene_Name)), "\n")
cat("Unique traits:", length(unique(plot_data$Trait)), "\n")

# Create the plot
cat("\n=== CREATING DOT PLOT ===\n")

p <- ggplot(plot_data, aes(x = Trait, y = Gene_Name)) +
  # Add points with size based on -log10(p-value) and color based on effect size
  geom_point(aes(
    size = log10_pvalue,
    color = Effect
  ), alpha = 0.8) +
  
  # Color scale: red for negative effects, green for positive effects
  scale_color_gradient2(
    low = "red", 
    mid = "white", 
    high = "darkgreen",
    midpoint = 0,
    name = "Effect Size",
    guide = guide_colorbar(title.position = "top")
  ) +
  
  # Size scale for -log10(p-value)
  scale_size_continuous(
    name = "-log10(p-value)",
    range = c(1, 8),  # Size range for dots
    guide = guide_legend(
      title.position = "top",
      override.aes = list(color = "black")
    )
  ) +
  
  # Customize theme
  theme_minimal() +
  theme(
    # Rotate x-axis labels
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    
    # Legend positioning
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    
    # Plot title
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    
    # Grid lines
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_line(color = "gray95", size = 0.3),
    
    # Plot margins
    plot.margin = margin(20, 20, 20, 20)
  ) +
  
  # Labels
  labs(
    title = "SNP-Trait Associations",
    subtitle = "Dot size: -log10(p-value), Color: Effect size",
    x = "Trait",
    y = "Gene",
    caption = paste("Showing", nrow(plot_data), "significant associations (p < 0.05)")
  )

# Display the plot
print(p)

# Save the plot
cat("\n=== SAVING PLOT ===\n")

# Save as high-resolution PNG
ggsave("SNP_Trait_Associations_Plot.png", 
       plot = p, 
       width = 12, 
       height = 8, 
       dpi = 300, 
       bg = "white")

# Save as PDF (vector format)
ggsave("SNP_Trait_Associations_Plot.pdf", 
       plot = p, 
       width = 12, 
       height = 8, 
       bg = "white")

# Save as SVG (for publications)
ggsave("SNP_Trait_Associations_Plot.svg", 
       plot = p, 
       width = 12, 
       height = 8, 
       bg = "white")

cat("Plots saved as:\n")
cat("- SNP_Trait_Associations_Plot.png (high-res)\n")
cat("- SNP_Trait_Associations_Plot.pdf (vector)\n")
cat("- SNP_Trait_Associations_Plot.svg (publication)\n")

# Create a summary table
cat("\n=== CREATING SUMMARY TABLE ===\n")

summary_table <- plot_data %>%
  arrange(P_Value) %>%
  select(Trait, Gene_Name, Effect, P_Value, log10_pvalue, Distance_To_Gene) %>%
  mutate(
    P_Value = formatC(P_Value, format = "e", digits = 2),
    Effect = round(Effect, 3),
    log10_pvalue = round(log10_pvalue, 2)
  )

# Save summary table
write.table(summary_table, "SNP_Trait_Associations_Summary.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Summary table saved as: SNP_Trait_Associations_Summary.txt\n")

# Show top associations
cat("\nTop 10 associations by significance:\n")
print(head(summary_table, 10))

# Show statistics
cat("\n=== PLOT STATISTICS ===\n")
cat("Associations plotted:", nrow(plot_data), "\n")
cat("Unique genes:", length(unique(plot_data$Gene_Name)), "\n")
cat("Unique traits:", length(unique(plot_data$Trait)), "\n")
cat("Effect size range:", round(min(plot_data$Effect), 3), "to", round(max(plot_data$Effect), 3), "\n")
cat("-log10(p-value) range:", round(min(plot_data$log10_pvalue), 2), "to", round(max(plot_data$log10_pvalue), 2), "\n")

cat("\n=== PLOT GENERATION COMPLETE ===\n")
cat("Your SNP-Trait associations plot is ready!\n")
