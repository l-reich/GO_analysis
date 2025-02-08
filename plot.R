library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)

# Define file path
data_path <- "plotting_data/"

# Define roots (namespaces)
roots <- c("molecular_function", "cellular_component", "biological_process")

# Define metrics
metrics <- c("gene_numbers", "gene_numbers_filt", "path_lengths", "size_diffs", "sig_genes_all", "sig_genes_filt")

# Create an empty list to store data
data_list <- list()

for (root in roots) {
  for (metric in metrics) {
    file_name <- paste0(data_path, root, "_", metric, ".txt")
    
    if (file.exists(file_name)) {
      # Read data as a numeric vector
      values <- scan(file_name, what = integer(), quiet = TRUE)
      
      # Store in a dataframe
      data_list[[paste(root, metric, sep = "_")]] <- data.frame(
        value = values,
        namespace = root,
        metric = metric
      )
    }
  }
}

# Combine all data frames into one
plot_data <- bind_rows(data_list)

# Define custom titles for each metric
custom_titles <- c(
  "gene_numbers" = "Gene Set Sizes",
  "gene_numbers_filt" = "Filtered Gene Set Sizes",
  "path_lengths" = "Path Lengths (Leaves to Root)",
  "size_diffs" = "Size Differences (Child-Parent)",
  "sig_genes_all" = "Significant Genes in Sets",
  "sig_genes_filt" = "Significant Genes in Filtered Sets"
)

ggplot(plot_data, aes(x = value, fill = namespace)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ metric, scales = "free", labeller = as_labeller(custom_titles)) +  # Custom facet titles
  scale_fill_manual(values = c("molecular_function" = "blue", "cellular_component" = "green", "biological_process" = "red")) +
  theme_minimal() +
  labs(x = "Counts", y = "Density") +  # No overall plot title
  theme(legend.position = "top")

#######
# plots from output
#######

# Define file paths
data_files <- list(
  "biological_process" = "comp_outs/biological_process_out.txt",
  "molecular_function" = "comp_outs/molecular_function_out.txt",
  "cellular_component" = "comp_outs/cellular_component_out.txt"
)

# Load and combine data
data_list <- lapply(names(data_files), function(namespace) {
  dt <- fread(data_files[[namespace]], sep = "\t")
  dt[, namespace := namespace]
})

data <- rbindlist(data_list)

# Convert numeric columns
data[, c("ks_stat", "ks_fdr", "size", "noverlap") := lapply(.SD, as.numeric), .SDcols = c("ks_stat", "ks_fdr", "size", "noverlap")]

# Save all plots and tables in a single PDF
pdf("gene_set_analysis.pdf", width = 12, height = 8)

## 1. Distribution of Enrichment Scores and P-values
ggplot(data, aes(x = ks_stat, fill = namespace)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ namespace, scales = "free") +
  theme_minimal() +
  labs(title = "Distribution of Enrichment Scores", x = "KS Statistic", y = "Density")

ggplot(data, aes(x = ks_fdr, fill = namespace)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ namespace, scales = "free") +
  theme_minimal() +
  labs(title = "Distribution of Adjusted P-values (FDR)", x = "KS FDR", y = "Density")

## 2. Scatter Plots
ggplot(data, aes(x = size, y = ks_stat, color = namespace)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ namespace, scales = "free") +
  theme_minimal() +
  labs(title = "KS Statistic vs. Gene Set Size", x = "Gene Set Size", y = "KS Statistic")

ggplot(data, aes(x = size, y = ks_fdr, color = namespace)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ namespace, scales = "free") +
  theme_minimal() +
  labs(title = "KS FDR vs. Gene Set Size", x = "Gene Set Size", y = "KS FDR")

## 3. Sorted Tables (Top 10 gene sets)
for (metric in c("noverlap", "ks_stat")) {
  for (ns in unique(data$namespace)) {
    cat("\n\nTop 5 Gene Sets by", metric, "for", ns, "\n")
    print(data[namespace == ns][order(get(metric), decreasing = TRUE)][1:5, .(term, name, size, noverlap, ks_stat, ks_fdr)])
  }
}

for (ns in unique(data$namespace)) {
  cat("\n\nTop 5 Gene Sets by adjusted p-value for", ns, "\n")
  print(data[namespace == ns][order(get("ks_fdr"), decreasing = FALSE)][1:5, .(term, name, size, noverlap, ks_stat, ks_fdr)])
}

dev.off()

cat("Analysis complete. Results saved in 'gene_set_analysis.pdf'\n")