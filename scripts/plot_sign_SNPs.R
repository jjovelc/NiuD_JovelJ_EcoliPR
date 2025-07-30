library(ggplot2)
library(tidyverse)

setwd("/Users/juanjovel/OneDrive/jj/UofC/data_analysis/dongyanNiu/BICOSS")

# 1. First, let's combine all results into one dataframe
result_files <- list.files(pattern="BICOSS_.*_significant_SNPs.txt")
all_results <- do.call(rbind, lapply(result_files, function(x) {
  df <- read.table(x, header=TRUE, sep="\t")
  if(nrow(df) > 0) return(df)
}))

# Convert p-values to -log10
all_results$neg_log_p <- -log10(all_results$P_Value)

# Create Manhattan-style plot
manhattan_plot <- ggplot(all_results, aes(x=SNP_Index, y=neg_log_p, color=Trait)) +
  geom_point(size=3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") +
  labs(x="SNP Position", y="-log10(P-value)", 
       title="Manhattan Plot of SNP Associations") +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="red")

# Effect size plot
effect_plot <- ggplot(all_results, aes(x=reorder(SNP_ID, Effect), y=Effect, fill=Trait)) +
  geom_bar(stat="identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8),
        legend.position = "right") +
  labs(x="SNP", y="Effect Size", 
       title="Effect Sizes of Significant SNPs")

# Count SNPs per trait
snp_counts <- all_results %>%
  group_by(Trait) %>%
  summarize(Count = n())

count_plot <- ggplot(snp_counts, aes(x=reorder(Trait, -Count), y=Count)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x="Trait", y="Number of Significant SNPs",
       title="Number of Significant SNPs per Trait")

# Save plots
ggsave("manhattan_plot.png", manhattan_plot, width=12, height=8)
ggsave("effect_plot.png", effect_plot, width=12, height=8)
ggsave("count_plot.png", count_plot, width=10, height=6)

# Create a heatmap of SNP-trait associations
# First, create a matrix of effects
effect_matrix <- reshape2::dcast(all_results, SNP_ID ~ Trait, value.var="Effect")
rownames(effect_matrix) <- effect_matrix$SNP_ID
effect_matrix$SNP_ID <- NULL

# Create heatmap
library(pheatmap)
png("association_heatmap.png", width=1200, height=800)
pheatmap(effect_matrix,
         display_numbers=TRUE,
         number_format="%.2f",
         main="Heatmap of SNP-Trait Associations",
         angle_col=45,
         fontsize=10)
dev.off()

# Add summary statistics
cat("\nSummary Statistics:\n")
cat("Total number of unique SNPs:", length(unique(all_results$SNP_ID)), "\n")
cat("Total number of traits with associations:", length(unique(all_results$Trait)), "\n")
cat("\nTop 5 strongest associations by p-value:\n")
print(head(all_results[order(all_results$P_Value),], 5))

# Add -log10(p-value)
all_results$neg_log_p <- -log10(all_results$P_Value)

# Create dot plot
# Create dot plot with larger axis labels
ggplot(all_results, aes(x=Trait, y=SNP_ID)) +
  geom_point(aes(size=neg_log_p, color=Effect)) +
  scale_size_continuous(name="-log10(p-value)", range=c(2, 10)) +
  scale_color_gradient2(low="red", mid="white", high="blue", 
                        midpoint=0, name="Effect Size") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=16),  # Doubled size
        axis.text.y=element_text(size=16),                      # Doubled size
        axis.title.x=element_text(size=20),                     # Doubled size
        axis.title.y=element_text(size=20),                     # Doubled size
        panel.grid.major=element_line(color="grey90"),
        legend.position="right",
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        plot.title=element_text(size=18)) +
  labs(x="Trait", y="SNP", 
       title="SNP-Trait Associations",
       subtitle="Dot size: -log10(p-value), Color: Effect size") +
  guides(size=guide_legend(order=1),
         color=guide_colorbar(order=2))

# Save the plot with higher resolution
ggsave("snp_trait_dotplot_red-blue.png", width=14, height=12, dpi=300)
