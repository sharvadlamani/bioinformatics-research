library(qqman)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Set the directory containing the CSV files
data_dir <- "/gpfs/commons/home/svadlamani/skat_all_races/AMR"
cell_directories <- c("coding", "LHX2", "NEUN", "OLIG2", "peripheralPU1nuclei")

# Initialize an empty dataframe to store combined results
combined_results <- data.frame()

# Loop through each cell and chromosome
for (dir in cell_directories) {
  for (chr_num in 1:22) {
    # Construct the file path for the current chromosome
    file_path <- file.path(data_dir, dir, "1", paste0("chr", chr_num, ".csv"))
    print(file_path)
    
    # Try to read the CSV file
    tryCatch({
      results <- read.csv(file_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      print(head(results))
      
      colnames(results) <- c("Gene", "P_Value1", "P_Value2", "P_Value3")
      results[, 2:4] <- lapply(results[, 2:4], function(x) as.numeric(as.character(x)))
      
      results$CHR <- chr_num
      results$Cell_Type <- dir
      
      results$BP <- sample(1:1e6, nrow(results), replace = TRUE)
      
      combined_results <- rbind(combined_results, results)
    }, error = function(e) {
      message(paste("File for chromosome", chr_num, "not found. Moving to the next chromosome."))
    })
  }
}

don <- combined_results %>%
  filter(!is.na(P_Value1), !is.na(P_Value2), !is.na(P_Value3)) %>%
  mutate(P1 = P_Value1, P2 = P_Value2, P3 = P_Value3) %>%
  select(CHR, BP, Gene, P1, P2, P3, Cell_Type)

if (nrow(don) == 0) {
  message("No valid data points remaining after filtering. Plot cannot be generated.")
} else {
  offset <- 5
  don <- don %>%
    group_by(CHR) %>%
    mutate(BPcum = BP + max(BP) * (CHR - 1) + offset * (CHR - 1)) %>%
    ungroup()
  
  axisdf <- don %>%
    group_by(CHR) %>%
    summarize(center = mean(BPcum))
  
  cell_colors <- c("red", "blue", "green", "orange", "purple")
  names(cell_colors) <- cell_directories
  
  # Set the annotation threshold
  annotation_threshold <- 2
  
  # Add annotation labels based on the threshold
  don <- don %>%
    mutate(Label = ifelse(-log10(P1) > annotation_threshold, as.character(Gene), NA))
  
  annotateTop <- TRUE  # Set to TRUE or FALSE as required
  
  if (!is.null(axisdf)) {
    plot <- ggplot(don, aes(x = BPcum, y = -log10(P1), color = Cell_Type)) +
      geom_point(alpha = 0.8, size = 1.3) +
      scale_color_manual(values = cell_colors) +
      scale_x_continuous(labels = axisdf$CHR, breaks = axisdf$center) +
      scale_y_continuous(expand = c(0, 0)) +
      xlab("Chromosome") +
      ylab("-log10(P-value)") +
      theme_bw() +
      theme(legend.position = "right",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            plot.margin = margin(20, 20, 50, 20, "pt"))  # Adjust top margin
    
    # Conditionally add annotations
    if (annotateTop) {
      plot <- plot +
        geom_text_repel(data = don %>% filter(!is.na(Label)), 
                        aes(label = Label), 
                        size = 3, 
                        box.padding = 0.1,    # Adjust padding
                        point.padding = 0.1,  # Adjust padding
                        max.iter = 200)       # Increase iterations for better positioning
    }
    
    ggsave("manhattan_plot_AMR.png", plot, dpi = 300, width = 20, height = 10)
    
    print(plot)
  }
  
  p_values_burden <- don$P1
  p_values_skat <- don$P2
  p_values_skato <- don$P3
  
  pdf("QQ_plots_AMR.pdf")
  
  qq(p_values_burden, main = "Burden QQ Plot AMR")
  qq(p_values_skat, main = "SKAT QQ Plot AMR")
  qq(p_values_skato, main = "SKATO QQ Plot AMR")
  
  dev.off()
}
