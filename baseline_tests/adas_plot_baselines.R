library(UpSetR)
library(qqman)
library(ggrepel)
library(dplyr)


# load in baseline results
results <- read.csv("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/data/baseline_gruyere.csv")
results$Type = results$Cell
results$Variant.Type = results$Cell


results$chr = as.numeric(gsub("[^0-9]", "", results$Chromosome))
color_column <- results$Type
#colors <- rainbow(length(unique(color_column)))
#colors <- c('purple', 'darkorange','red', 'darkgreen', 'blue4')
#colors <- c('#984EA3', '#FF7F00', '#E41A1C',  '#4DAF4A', '#377EB8')
colors <- c('#8491B4FF', '#F39B7FFF', '#DC0000FF',  '#00A087FF', '#3C5488FF')
colors <- c('#69799EFF', '#D37247FF', '#DC0000FF',  '#00A087FF', '#3C5488FF')



color_vector <- colors[as.numeric(factor(color_column))]
results$color = color_vector

TESTS = colnames(results)[c(3,4,5,7,8,9,10,11,12,13,14)]
THRESHOLD = 0.00001

for(TEST in c("gruyere")){
  df = results[complete.cases(results[[TEST]]), c('Gene','Variant.Type','chr','bp', TEST, 'color')]
  colnames(df) <- c('SNP','VARIANT.TYPE','CHR','BP', 'P', 'COLOR')
  significant = df[df[[TEST]]<=THRESHOLD,]
  gwasResults = df
  
  gwasResults %>% 
    filter(-log10(P)>1)
  
  don <- gwasResults %>% 
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>%
    
    # Add highlight and annotation information
    mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
    mutate( is_annotate=ifelse(-log10(P)>-log10(THRESHOLD), "yes", "no")) 
  
  # Prepare X axis
  axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  # Make the plot
  plot = ggplot(don, aes(x=BPcum, y=-log10(P))) +
    
    # Show all points
    #geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3, show.legend = FALSE) +
    #scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    geom_point( aes(color=VARIANT.TYPE), alpha=0.8, size=1.3, show.legend = TRUE) +
    scale_color_manual(values = rep(c(colors), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP),colour = subset(don, is_annotate=="yes")$COLOR, size=3, show.legend=FALSE, max.overlaps=100) +
    xlab("Chromosome") + ylim(0, 11)+
    #geom_point(aes(color = VARIANT.TYPE)) +  theme(legend_text=element_text(size=14)) + 
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="right",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
   
  
  #ggsave(paste0("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/figures/baseline_figs/", TEST, '_manhattan_c.png'), plot, dpi = 600, width = 13, height = 8.14)
}

print(plot)
ggsave("/gpfs/commons/groups/knowles_lab/data/ADSP_reguloML/Alzheimer-RV/figures/baseline_figs/gruyere_manhattan.png", plot, dpi = 600, width = 13, height = 8)

#manhattan(df, chr = "CHR", bp='BP', snp='SNP', p = 'P', annotatePval = THRESHOLD, col = c('grey','lightblue'))

#ggplot(results, aes(x =Inf, y = Inf, color = Variant.Type), show.legend = TRUE) + geom_point() + scale_color_manual(values = colors) + theme(legend_text=element_text(size=14))







