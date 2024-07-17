# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ggrepel)

# Set the directories containing the CSV files for AFR and AMR
data_dir_afr <- "/gpfs/commons/home/svadlamani/skat_all_maf_range/0.001"
data_dir_amr <- "/gpfs/commons/home/svadlamani/skat_all_maf_range/5e-05"
cell_directories <- c("coding", "LHX2", "NEUN", "OLIG2", "peripheralPU1nuclei")

# Function to read data from CSV files
read_data <- function(data_dir, cell_directories, label) {
  combined_results <- data.frame()
  
  for (dir in cell_directories) {
    for (chr_num in 1:22) {
      file_path <- file.path(data_dir, dir, "1", paste0("chr", chr_num, ".csv"))
      print(file_path)
      
      tryCatch({
        results <- read.csv(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        print(head(results))
        
        colnames(results) <- c("Gene", "P_Value1", "P_Value2", "P_Value3")
        results[, 2:4] <- lapply(results[, 2:4], function(x) as.numeric(as.character(x)))
        
        results$CHR <- chr_num
        results$Cell_Type <- dir
        results$Label <- label
        
        results$BP <- sample(1:1e6, nrow(results), replace = TRUE)
        
        combined_results <- rbind(combined_results, results)
      }, error = function(e) {
        message(paste("File for chromosome", chr_num, "not found. Moving to the next chromosome."))
      })
    }
  }
  
  return(combined_results)
}

# Read data from both AFR and AMR directories
combined_results_afr <- read_data(data_dir_afr, cell_directories, "AFR")
combined_results_amr <- read_data(data_dir_amr, cell_directories, "AMR")

# Filter and prepare the data for plotting
prepare_data <- function(combined_results) {
  don <- combined_results %>%
    filter(!is.na(P_Value1) & !is.na(P_Value2) & !is.na(P_Value3)) %>%
    mutate(P1 = -log10(P_Value1), P2 = -log10(P_Value2), P3 = -log10(P_Value3)) %>%
    group_by(Gene, Cell_Type) %>%
    summarize(P1 = min(P1), P2 = min(P2), P3 = min(P3), .groups = 'drop') # Ensure unique data points for each gene and cell type
  
  return(don)
}

don_afr <- prepare_data(combined_results_afr)
don_amr <- prepare_data(combined_results_amr)

# Merge data and prepare for plotting
merged_data <- merge(don_afr, don_amr, by = c("Gene", "Cell_Type"), suffixes = c("_AFR", "_AMR"))

# Debugging: Check the merged data
head(merged_data)
sum(duplicated(merged_data))

# Function to get top genes based on P-value threshold
get_top_genes <- function(df, threshold, column) {
  df %>%
    arrange(desc(.data[[column]])) %>%
    filter(.data[[column]] > 4) %>%  # Filter genes with -log(p value) > 4
    head(threshold) %>%
    mutate(Label = as.character(Gene))
}

# Function to plot scatter plot with outlier labels and color coding
plot_scatter <- function(merged_data, test, test_label) {
  p_val_afr <- paste0(test, "_AFR")
  p_val_amr <- paste0(test, "_AMR")
  
  # Filter out rows with NA, NaN, or Inf values in both columns
  merged_data <- merged_data %>%
    filter(!is.na(.data[[p_val_afr]]) & !is.nan(.data[[p_val_afr]]) & !is.infinite(.data[[p_val_afr]]),
           !is.na(.data[[p_val_amr]]) & !is.nan(.data[[p_val_amr]]) & !is.infinite(.data[[p_val_amr]]))
  
  top_genes_afr <- get_top_genes(merged_data, 50, p_val_afr)
  top_genes_amr <- get_top_genes(merged_data, 50, p_val_amr)
  
  top_genes <- unique(c(top_genes_afr$Gene, top_genes_amr$Gene))
  
  merged_data$Label <- ifelse(merged_data$Gene %in% top_genes | 
                                merged_data[[p_val_afr]] >= 2 |
                                merged_data[[p_val_amr]] >= 2, 
                              as.character(merged_data$Gene), NA)
  
  # Calculate R-squared
  lm_model <- lm(merged_data[[p_val_afr]] ~ merged_data[[p_val_amr]], data = merged_data)
  r_squared <- summary(lm_model)$r.squared
  
  # Append R-squared to the plot title
  plot_title <- paste0("Scatter plot of -log10(P-values) for ", test_label, "\n", 
                       "R-squared: ", round(r_squared, digits = 4))
  
  # Create the scatter plot
  plot <- ggplot(merged_data, aes(x = .data[[p_val_amr]], y = .data[[p_val_afr]], color = Cell_Type)) +
    geom_point(alpha = 0.6) +
    geom_text_repel(data = merged_data %>% filter(!is.na(Label)), aes(label = Label), size = 3, box.padding = 0.3,
                    max.overlaps = 28) +  # Adjust max.overlaps to allow more labels
    labs(x = paste0("-log10(P-value) ", test_label, " 5e-05"),
         y = paste0("-log10(P-value) ", test_label, " 0.001"),
         title = plot_title,
         color = "Cell Type") +
    theme_minimal()
  
  print(plot)
  
  ggsave(paste0("scatter_plot_", tolower(test_label), ".png"), plot = plot, dpi = 300, width = 10, height = 8)
}

# Plot scatter plots for Burden (P1), SKAT (P2), and SKATO (P3)
plot_scatter(merged_data, "P1", "Burden_0.001_5e-05")
plot_scatter(merged_data, "P2", "SKAT_0.001_5e-05")
plot_scatter(merged_data, "P3", "SKATO_0.001_5e-05")