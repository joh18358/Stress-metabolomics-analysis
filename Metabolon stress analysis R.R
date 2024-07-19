rm(list=ls())

library(gplots)
library(ggplot2)
library(VennDiagram)
library(viridis)
library(readxl)
library(conflicted)
library(dplyr)
library(tidyverse) # contains ggplot for volcano plots
library(RColorBrewer) # to color plots
library(ggrepel) # annotating volcano plot
library(ggfortify) # PCA plots
library(factoextra) # PCA plot ellipses
library(MSnSet.utils) # Different PCA plot option (simpler?)
library('corrr') # PCA plot 3rd option
library(ggcorrplot) # PCA plot 3rd option
library("FactoMineR") # PCA plot 3rd option
library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
library("xlsx")
library(pheatmap)

if(!require('DiscriMiner')) {
  install.packages('DiscriMiner')
  library('DiscriMiner')
}

conflicts_prefer(base::setdiff, base::intersect, dplyr::filter())


setwd("C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel")

# Upload data files ####

#PCA file
PCA_file <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/PCA_file_metabolon_stress.xlsx"))

# Metabolite lists from 2-group comparisons
mdxvsWT_baseline_sig_metabolites <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/mdxvsWT_baseline_DCMs.xlsx"))

mdxvsMTBD_baseline_sig_metabolites <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/mdxvsMTBD_baseline_DCMs.xlsx"))

mdxvsWT_scruff_sig_metabolites <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/mdxvsWT_scruff_DCMs.xlsx"))

mdxvsMTBD_scruff_sig_metabolites <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/mdxvsMTBD_scruff_DCMs.xlsx"))

mdx_scruffvsbaseline_sig_metabolites <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/mdx_scruffvsbaseline_DCMs.xlsx"))

WT_scruffvsbaseline_sig_metabolites <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/WT_scruffvsbaseline_DCMs.xlsx"))

MTBD_scruffvsbaseline_sig_metabolites <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/MTBD_scruffvsbaseline_DCMs.xlsx"))

MTBDvsWT_baseline_sig_metabolites <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/MTBDvsWT_baseline_DCMs.xlsx"))

MTBDvsWT_scruff_sig_metabolites <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/MTBDvsWT_scruff_DCMs.xlsx"))

# Heatmap Data

heatmap_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/heatmap_avg_data_file.xlsx"))

# Volcano plot data

mdxvsWT_baseline_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/mdxvsWT_baseline_volcano.xlsx"))

mdxvsMTBD_baseline_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/mdxvsMTBD_baseline_volcano.xlsx"))

mdxvsWT_scruff_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/mdxvsWT_scruff_volcano.xlsx"))

mdxvsMTBD_scruff_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/mdxvsMTBD_scruff_volcano.xlsx"))

mdx_scruffvsbaseline_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/mdx_scruffvsbaseline_volcano.xlsx"))

WT_scruffvsbaseline_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/WT_scruffvsbaseline_volcano.xlsx"))

MTBD_scruffvsbaseline_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/MTBD_scruffvsbaseline_volcano.xlsx"))

MTBDvsWT_baseline_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/MTBDvsWT_baseline_volcano.xlsx"))

MTBDvsWT_scruff_volcano_data <-
  data.frame(read_excel(
    "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/R files metabolon/MTBDvsWT_scruff_volcano.xlsx"))

----

  
####### PCA plot ####### #############

PCA_no_legend <- PCA_file[2:73]

PCA_no_missing <- PCA_no_legend[ , colSums(is.na(PCA_no_legend))==0]

pca_res_metabolon <- prcomp(PCA_no_missing, scale. = TRUE)

PCi_metabolon<-data.frame(pca_res_metabolon$x,Legend=PCA_file$Legend)

ggplot(PCi_metabolon,aes(x=PC1,y=PC2,col=Legend,frame=T))+
  geom_point(size=3,alpha=0.6) + #Size and alpha just for fun
  theme_minimal() +
  theme(text = element_text(size=15)) +
  stat_ellipse(alpha=0.03,geom = "polygon", aes(fill=Legend)) +
  scale_color_manual(values = c("#414487","#7ad151","#2a788e","#fde725","#440154","#22a884"))

## Other PCA plot option using corrr and FactoMineR ##

corr_matrix_metabolon <- cor(PCA_no_missing)
#ggcorrplot(corr_matrix_ECF)

data_normalized_metabolon <- scale(PCA_no_missing)
head(data_normalized_metabolon)

# For some reason, can't use princomp like example showed
# Something about needing more units than variables and R-mode vs Q-mode
data.pca_metabolon <- prcomp(data_normalized_metabolon)
summary(data.pca_metabolon)

# Scree plot
fviz_eig(data.pca_metabolon, addlabels = TRUE)

# Biplot
fviz_pca_var(data.pca_metabolon, col.var = "cos2",
             gradient.cols=c("black","turquoise","purple"),
             ggtheme = theme_minimal(),
             select.var = list(name = NULL, cos2 = NULL, contrib=5),
             repel=TRUE)

cos2_metabolon <- fviz_cos2(data.pca_metabolon, choice = "var", axes = 1:2, 
                        sort.val = "desc",
                        top=20) # Can change top to any numerical value
cos2_metabolon_PC1only <- fviz_cos2(data.pca_metabolon, choice = "var", axes = 1, 
                                sort.val = "desc",
                                top=20) # Can change top to any numerical value

topfeatures_metabolon <- as.data.frame(cos2_metabolon$data)
topfeatures_metabolon_PC1 <- as.data.frame(cos2_metabolon_PC1only$data)

#write.xlsx(topfeatures_metabolon, 
#           "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/metabolon_PC1_PC2_top_features.xlsx")
#write.xlsx(topfeatures_metabolon_PC1, 
#          "C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/metabolon_PC1_top_features.xlsx")

# Most useful and visually appealing biplot
# Can change contrib number to add or subtract vectors
fviz_pca_biplot(data.pca_metabolon, axes = c(1,2), label = "var", repel = TRUE,
                col.var = "black",
                geom = "point",
                pointsize = 2,
                habillage=PCA_file$Legend,
                palette = c("#414487","#7ad151","#2a788e","#fde725","#440154","#22a884"), 
                addEllipses=TRUE, ellipse.level=0.95,
                ellipse.type = "t",
                select.var = list(contrib = 15),
                ggtheme = theme_minimal())
---
  
  
  

##### 4-way Venn Diagram Metabolon data
###############
# Venn Diagrams

#library(RColorBrewer)
#myCol <- brewer.pal(3, "Pastel2")

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
# Suppresses data file generation each time a Venn diagram is generated
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

#####
Metabolon_4D_venn <- display_venn(Metabolon_4D_venn_list <- list(
  'mdx vs WT baseline' = mdxvsWT_baseline_sig_metabolites$'Metabolite',
  'mdx vs MTBD baseline' = mdxvsMTBD_baseline_sig_metabolites$'Metabolite',
  'mdx vs WT post-scruff' = mdxvsWT_scruff_sig_metabolites$'Metabolite',
  'mdx vs MTBD post-scruff' = mdxvsMTBD_scruff_sig_metabolites$'Metabolite'), 
  # Circles
  lwd = 2,
  lty = 'blank',
  category.names = c(
    "mdx vs WT baseline DCM",
    "mdx vs MTBD baseline DCM",
    "mdx vs WT post-scruff DCM", 
    "mdx vs MTBD post-scruff DCM"),
  fill = c("#776D5A","#987D7C","#A09CB0", "#A3B9C9"),
  # Numbers
  cex = 1.2,
  fontface = "italic",
  # Set names
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-9, 2, -20, -20),
  cat.dist = c(0.135, .13, 0.055, 0.05),
  main = "Circulating Stress Metabolite Comparisons WT, mdx, and MTBD",
  main.fontface = "bold"
) 

-----

###### Group overlap/non-overlap comparisons #######
# Baseline non-overlap 2-group DCMs #########

mdxvsWT_baseline_not_MTBD_baseline <- 
  data.frame('mdxvsWT_baseline_not_MTBD_baseline' = setdiff
             (mdxvsWT_baseline_sig_metabolites$'Metabolite',
               mdxvsMTBD_baseline_sig_metabolites$'Metabolite')) %>% 
  'colnames<-' ("Metabolite")

mdxvsWT_baseline_unique <- 
  data.frame('mdxvsWT_baseline_unique' = setdiff
             (mdxvsWT_baseline_not_MTBD_baseline$'Metabolite',
               mdxvsWT_scruff_sig_metabolites$'Metabolite')) %>% 
  'colnames<-' ("Metabolite")

# write.xlsx(mdxvsWT_baseline_unique, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/mdxvsWT_baseline_only_DCMs.xlsx")

mdxvsMTBD_baseline_not_WT_baseline <- 
  data.frame('mdxvsMTBD_baseline_not_WT_baseline' = setdiff
             (mdxvsMTBD_baseline_sig_metabolites$'Metabolite',
               mdxvsWT_baseline_sig_metabolites$'Metabolite')) %>% 
  'colnames<-' ("Metabolite")

mdxvsMTBD_baseline_unique <- 
  data.frame('mdxvsMTBD_baseline_unique' = setdiff
             (mdxvsMTBD_baseline_not_WT_baseline$'Metabolite',
               mdxvsWT_scruff_sig_metabolites$'Metabolite')) %>% 
  'colnames<-' ("Metabolite")

# write.xlsx(mdxvsMTBD_baseline_unique, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/mdxvsMTBD_baseline_only_DCMs.xlsx")
-----

# Scruff non-overlap 2-group DCMs #######
mdxvsWT_scruff_not_MTBD_scruff <- 
  data.frame('mdxvsWT_scruff_not_MTBD_scruff' = setdiff
             (mdxvsWT_scruff_sig_metabolites$'Metabolite',
               mdxvsMTBD_scruff_sig_metabolites$'Metabolite')) %>% 
  'colnames<-' ("Metabolite")

mdxvsWT_scruff_unique <- 
  data.frame('mdxvsWT_scruff_unique' = setdiff
             (mdxvsWT_scruff_not_MTBD_scruff$'Metabolite',
               mdxvsWT_baseline_sig_metabolites$'Metabolite')) %>% 
  'colnames<-' ("Metabolite")

# write.xlsx(mdxvsWT_scruff_unique, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/mdxvsWT_scruff_only_DCMs.xlsx")

mdxvsMTBD_scruff_not_WT_scruff <- 
  data.frame('mdxvsMTBD_scruff_not_WT_scruff' = setdiff
             (mdxvsMTBD_scruff_sig_metabolites$'Metabolite',
               mdxvsWT_scruff_sig_metabolites$'Metabolite')) %>% 
  'colnames<-' ("Metabolite")

mdxvsMTBD_scruff_unique <- 
  data.frame('mdxvsMTBD_scruff_unique' = setdiff
             (mdxvsMTBD_scruff_not_WT_scruff$'Metabolite',
               mdxvsMTBD_baseline_sig_metabolites$'Metabolite')) %>% 
  'colnames<-' ("Metabolite")

# write.xlsx(mdxvsMTBD_scruff_unique, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/mdxvsmTBD_scruff_only_DCMs.xlsx")

-----

# Baseline overlap DCMS mdx vs. WT and mdx vs. MTBD #####

baseline_overlap_DCM <- 
  data.frame('baseline_overlap_DCM' = intersect
             (mdxvsWT_baseline_sig_metabolites$'Metabolite',
               mdxvsMTBD_baseline_sig_metabolites$'Metabolite')) %>% 
  'colnames<-' ("Metabolite")

uniquely_overlapping_baseline_DCM <-
  baseline_overlap_DCM %>% 
  base::subset(!(.$'Metabolite' %in% mdxvsWT_scruff_sig_metabolites$'Metabolite')) %>%
  as.data.frame(.) %>%
  base::subset(!(.$'Metabolite' %in% mdxvsMTBD_scruff_sig_metabolites$'Metabolite')) %>%
  as.data.frame(.)

# write.xlsx(uniquely_overlapping_baseline_DCM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/overlapping_baseline_DCMs.xlsx")

-----

# Post-scruff overlap DCMS mdx vs. WT and mdx vs. MTBD ####

scruff_overlap_DCM <- 
  data.frame('scrff_overlap_DCM' = intersect
             (mdxvsWT_scruff_sig_metabolites$'Metabolite',
               mdxvsMTBD_scruff_sig_metabolites$'Metabolite')) %>% 
  'colnames<-' ("Metabolite")

uniquely_overlapping_scruff_DCM <-
  scruff_overlap_DCM %>% 
  base::subset(!(.$'Metabolite' %in% mdxvsWT_baseline_sig_metabolites$'Metabolite')) %>%
  as.data.frame(.) %>%
  base::subset(!(.$'Metabolite' %in% mdxvsMTBD_baseline_sig_metabolites$'Metabolite')) %>%
  as.data.frame(.)

# write.xlsx(uniquely_overlapping_scruff_DCM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/overlapping_scruff_DCMs.xlsx")

----

# Overlapping all groups #####
all_overlap_DCM <-
  data.frame('all_overlap_DCM' = intersect
      (scruff_overlap_DCM$'Metabolite',
          baseline_overlap_DCM$'Metabolite')) %>%
  'colnames<-' ("Metabolite")

#write.xlsx(all_overlap_DCM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/overlapping_baseline_and_scruff_DCMs.xlsx")

-----
  
# Unique Post-scruff vs baseline mdx DCMs ####

mdx_scruff_unique_DCM <-
  mdx_scruffvsbaseline_sig_metabolites %>% 
  base::subset(!(.$'Metabolite' %in% WT_scruffvsbaseline_sig_metabolites$'Metabolite')) %>%
  as.data.frame(.) %>%
  base::subset(!(.$'Metabolite' %in% MTBD_scruffvsbaseline_sig_metabolites$'Metabolite')) %>%
  as.data.frame(.)

mdx_scruff_unique_DCM <-
  mdx_scruff_unique_DCM %>% 
  base::subset(!(.$'Metabolite' %in% mdxvsWT_baseline_sig_metabolites$'Metabolite')) %>%
  as.data.frame(.) %>%
  base::subset(!(.$'Metabolite' %in% mdxvsMTBD_baseline_sig_metabolites$'Metabolite')) %>%
  as.data.frame(.)

# write.xlsx(mdx_scruff_unique_DCM, file="C:/Users/JOH18358/OneDrive/Ervasti Lab/Metabolon stress panel/unique_mdx_scruff_DCMs.xlsx")

-----
  
  
  

# Heat map ######

scaled_heatmap_data <- scale(heatmap_data[2:48])
row.names(scaled_heatmap_data) <- heatmap_data[,1]

unscaled_heatmap_data <- heatmap_data[2:48]
row.names(unscaled_heatmap_data) <- heatmap_data[,1]

clustered_heatmap_data <- pheatmap(scaled_heatmap_data, main="Heatmap")

row_clustered_heatmap_data <- pheatmap(unscaled_heatmap_data, cluster_cols = FALSE, 
                                       scale = 'column',
                                       main="Row Cluster Heatmap")

col_clustered_heatmap_data <- pheatmap(scaled_heatmap_data, cluster_rows = FALSE, 
                                       main="Column Cluster Heatmap")

not_clustered_heatmap_data <- pheatmap(scaled_heatmap_data, cluster_rows = FALSE, 
                                       cluster_cols = FALSE,
                                       cellheight = 8,
                                       fontsize_row = 9,
                                       #scale = 'column',
                                       main="No Clustering Heatmap")

matrix_heatmap_data <- as.matrix(scaled_heatmap_data)

heatmap(matrix_heatmap_data, Colv = NA, Rowv = NA, scale = "column")  

# Volcano Plots ######

# Defining Volcano plot function 
volcano_fun <- function(volcano_df){
  volcano_df %>%
    # Add a threshold for significant observations
    mutate(threshold =
             if_else(`L2FC` >= 0 & `log10pval` >= 1.3, "A",
                     if_else(`L2FC` <= 0 & `log10pval` >= 1.3,"B","C"))) %>%
    # Plot with points coloured according to the threshold
    ggplot(aes(x=`L2FC`, y=`log10pval`, colour = threshold)) +
    geom_point(alpha = 1, size = 3, shape = 16) + # Alpha sets the transparency of the points
    # Add dotted lines to indicate the threshold, semi-transparent
    geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) + 
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
    coord_cartesian(ylim = c(0, 5), xlim = c(-4, 5)) +
    labs(color = 'Legend',
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    # Set the colour of the points
    scale_colour_manual(values = c("A"= "#674961", "B"= "#80ADA0", 
                                   "C" = "#A2A7A5"),
                        labels = c("Upregulated", "Downregulated", "Not Significant")) +
    #xlab("log2(fold change)") + ylab("-log10(adjusted P-value)") + # Relabel the axes
    theme_minimal() + #Set the theme
    theme(legend.position = c(.98, .98),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6,),
          legend.box.background = element_rect(color="#427681",linewidth =1)) +
    theme(text = element_text(size=14)) 
}

# Volcano plot function if you want labels
volcano_fun_label <- function(volcano_df){
  volcano_df %>%
    # Add a threshold for significant observations
    mutate(threshold =
             if_else(`L2FC` >= 0 & `log10pval` >= 1.3, "A",
                     if_else(`L2FC` <= 0 & `log10pval` >= 1.3,"B","C"))) %>%
    # Plot with points coloured according to the threshold
    ggplot(aes(x=`L2FC`, y=`log10pval`, colour = threshold, label = `DCMlabel`)) +
    geom_point(alpha = 1, size = 3, shape = 16) + # Alpha sets the transparency of the points
    # Add dotted lines to indicate the threshold, semi-transparent
    geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) + 
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
    coord_cartesian(ylim = c(0, 5), xlim = c(-4, 5)) +
    labs(color = 'Legend',
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    # Set the colour of the points
    scale_colour_manual(values = c("A"= "#674961", "B"= "#80ADA0", 
                                   "C" = "#A2A7A5"),
                        labels = c("Upregulated", "Downregulated", "Not Significant")) +
    #xlab("log2(fold change)") + ylab("-log10(adjusted P-value)") + # Relabel the axes
    theme_minimal() + #Set the theme
    theme(legend.position = c(0.98, 0.98),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6,),
          legend.box.background = element_rect(color="#427681",linewidth =1)) +
    theme(text = element_text(size=14)) +
    geom_text_repel(max.overlaps = Inf, key_glyph = draw_key_point)
}

volcano_fun2 <- function(volcano_df){
  volcano_df %>%
    # Add a threshold for significant observations
    mutate(threshold =
             if_else(`L2FC` >= 0 & `log10pval` >= 1.3, "A",
                     if_else(`L2FC` <= 0 & `log10pval` >= 1.3,"B","C"))) %>%
    # Plot with points coloured according to the threshold
    ggplot(aes(x=`L2FC`, y=`log10pval`, colour = threshold)) +
    geom_point(alpha = 1, size = 3, shape = 16) + # Alpha sets the transparency of the points
    # Add dotted lines to indicate the threshold, semi-transparent
    geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) + 
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
    coord_cartesian(ylim = c(0, 8), xlim = c(-4, 5)) +
    labs(color = 'Legend',
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    # Set the colour of the points
    scale_colour_manual(values = c("A"= "#674961", "B"= "#80ADA0", 
                                   "C" = "#A2A7A5"),
                        labels = c("Upregulated", "Downregulated", "Not Significant")) +
    #xlab("log2(fold change)") + ylab("-log10(adjusted P-value)") + # Relabel the axes
    theme_minimal() + #Set the theme
    theme(legend.position = c(0.03, .98),
          legend.justification = c("left", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6,),
          legend.box.background = element_rect(color="#427681",linewidth =1)) +
    theme(text = element_text(size=14)) 
}
volcano_fun_label2 <- function(volcano_df){
  volcano_df %>%
    # Add a threshold for significant observations
    mutate(threshold =
             if_else(`L2FC` >= 0 & `log10pval` >= 1.3, "A",
                     if_else(`L2FC` <= 0 & `log10pval` >= 1.3,"B","C"))) %>%
    # Plot with points coloured according to the threshold
    ggplot(aes(x=`L2FC`, y=`log10pval`, colour = threshold, label = `DCMlabel`)) +
    geom_point(alpha = 1, size = 3, shape = 16) + # Alpha sets the transparency of the points
    # Add dotted lines to indicate the threshold, semi-transparent
    geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) + 
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.5) +
    coord_cartesian(ylim = c(0, 8), xlim = c(-4, 5)) +
    labs(color = 'Legend',
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
    # Set the colour of the points
    scale_colour_manual(values = c("A"= "#674961", "B"= "#80ADA0", 
                                   "C" = "#A2A7A5"),
                        labels = c("Upregulated", "Downregulated", "Not Significant")) +
    #xlab("log2(fold change)") + ylab("-log10(adjusted P-value)") + # Relabel the axes
    theme_minimal() + #Set the theme
    theme(legend.position = c(0.03, 0.98),
          legend.justification = c("left", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6,),
          legend.box.background = element_rect(color="#427681",linewidth =1)) +
    theme(text = element_text(size=14)) +
    geom_text_repel(max.overlaps = Inf, key_glyph = draw_key_point)
}

## volcano_data plot for mdx vs WT baseline

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
mdxvsWT_baseline_volcano_data$diffexpressed <- NA
# if log2Foldchange > 0 and -log10pval > 1.3, set as "UP"
mdxvsWT_baseline_volcano_data$diffexpressed[mdxvsWT_baseline_volcano_data$L2FC > 0 & 
                                         mdxvsWT_baseline_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < 0 and -log10pval > 1.3, set as "DOWN"
mdxvsWT_baseline_volcano_data$diffexpressed[mdxvsWT_baseline_volcano_data$L2FC < 0 & 
                                             mdxvsWT_baseline_volcano_data$log10pval > 1.3] <- "DOWN"

# Add column used to label points
mdxvsWT_baseline_volcano_data$DCMlabel <- ifelse(grepl("UP", mdxvsWT_baseline_volcano_data$diffexpressed), 
        mdxvsWT_baseline_volcano_data$Metabolite, 
        ifelse(grepl("DOWN", mdxvsWT_baseline_volcano_data$diffexpressed), 
               mdxvsWT_baseline_volcano_data$Metabolite, NA))

volcano_fun(mdxvsWT_baseline_volcano_data)
volcano_fun_label(mdxvsWT_baseline_volcano_data)

## volcano_data plot for mdx vs MTBD baseline

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
mdxvsMTBD_baseline_volcano_data$diffexpressed <- "NO"
# if log2Foldchange > 0 and -log10pval > 1.3, set as "UP"
mdxvsMTBD_baseline_volcano_data$diffexpressed[mdxvsMTBD_baseline_volcano_data$L2FC > 0 & 
                                         mdxvsMTBD_baseline_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < 0 and -log10pval > 1.3, set as "DOWN"
mdxvsMTBD_baseline_volcano_data$diffexpressed[mdxvsMTBD_baseline_volcano_data$L2FC < 0 & 
                                         mdxvsMTBD_baseline_volcano_data$log10pval > 1.3] <- "DOWN"

# Add column used to label points
mdxvsMTBD_baseline_volcano_data$DCMlabel <- ifelse(grepl("UP", mdxvsMTBD_baseline_volcano_data$diffexpressed), 
                                                 mdxvsMTBD_baseline_volcano_data$Metabolite, 
                                                 ifelse(grepl("DOWN", mdxvsMTBD_baseline_volcano_data$diffexpressed), 
                                                        mdxvsMTBD_baseline_volcano_data$Metabolite, NA))

volcano_fun(mdxvsMTBD_baseline_volcano_data)
volcano_fun_label(mdxvsMTBD_baseline_volcano_data)

## volcano_data plot for WT vs MTBD baseline

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
MTBDvsWT_baseline_volcano_data$diffexpressed <- "NO"
# if log2Foldchange > 0 and -log10pval > 1.3, set as "UP"
MTBDvsWT_baseline_volcano_data$diffexpressed[MTBDvsWT_baseline_volcano_data$L2FC > 0 & 
                                               MTBDvsWT_baseline_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < 0 and -log10pval > 1.3, set as "DOWN"
MTBDvsWT_baseline_volcano_data$diffexpressed[MTBDvsWT_baseline_volcano_data$L2FC < 0 & 
                                               MTBDvsWT_baseline_volcano_data$log10pval > 1.3] <- "DOWN"

# Add column used to label points
MTBDvsWT_baseline_volcano_data$DCMlabel <- ifelse(grepl("UP", MTBDvsWT_baseline_volcano_data$diffexpressed), 
                                                   MTBDvsWT_baseline_volcano_data$Metabolite, 
                                                   ifelse(grepl("DOWN", MTBDvsWT_baseline_volcano_data$diffexpressed), 
                                                          MTBDvsWT_baseline_volcano_data$Metabolite, NA))

volcano_fun(MTBDvsWT_baseline_volcano_data)
volcano_fun_label(MTBDvsWT_baseline_volcano_data)

## volcano_data plot for mdx vs WT post-scruff

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
mdxvsWT_scruff_volcano_data$diffexpressed <- NA
# if log2Foldchange > 0 and -log10pval > 1.3, set as "UP"
mdxvsWT_scruff_volcano_data$diffexpressed[mdxvsWT_scruff_volcano_data$L2FC > 0 & 
                                              mdxvsWT_scruff_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < 0 and -log10pval > 1.3, set as "DOWN"
mdxvsWT_scruff_volcano_data$diffexpressed[mdxvsWT_scruff_volcano_data$L2FC < 0 & 
                                              mdxvsWT_scruff_volcano_data$log10pval > 1.3] <- "DOWN"

# Add column used to label points
mdxvsWT_scruff_volcano_data$DCMlabel <- ifelse(grepl("UP", mdxvsWT_scruff_volcano_data$diffexpressed), 
                                                 mdxvsWT_scruff_volcano_data$Metabolite, 
                                                 ifelse(grepl("DOWN", mdxvsWT_scruff_volcano_data$diffexpressed), 
                                                        mdxvsWT_scruff_volcano_data$Metabolite, NA))

volcano_fun2(mdxvsWT_scruff_volcano_data)
volcano_fun_label2(mdxvsWT_scruff_volcano_data)

## volcano_data plot for mdx vs MTBD post-scruff

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
mdxvsMTBD_scruff_volcano_data$diffexpressed <- "NO"
# if log2Foldchange > 0 and -log10pval > 1.3, set as "UP"
mdxvsMTBD_scruff_volcano_data$diffexpressed[mdxvsMTBD_scruff_volcano_data$L2FC > 0 & 
                                                mdxvsMTBD_scruff_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < 0 and -log10pval > 1.3, set as "DOWN"
mdxvsMTBD_scruff_volcano_data$diffexpressed[mdxvsMTBD_scruff_volcano_data$L2FC < 0 & 
                                                mdxvsMTBD_scruff_volcano_data$log10pval > 1.3] <- "DOWN"

# Add column used to label points
mdxvsMTBD_scruff_volcano_data$DCMlabel <- ifelse(grepl("UP", mdxvsMTBD_scruff_volcano_data$diffexpressed), 
                                                   mdxvsMTBD_scruff_volcano_data$Metabolite, 
                                                   ifelse(grepl("DOWN", mdxvsMTBD_scruff_volcano_data$diffexpressed), 
                                                          mdxvsMTBD_scruff_volcano_data$Metabolite, NA))

volcano_fun(mdxvsMTBD_scruff_volcano_data)
volcano_fun_label(mdxvsMTBD_scruff_volcano_data)

## volcano_data plot for WT vs MTBD post-scruff

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
MTBDvsWT_scruff_volcano_data$diffexpressed <- "NO"
# if log2Foldchange > 0 and -log10pval > 1.3, set as "UP"
MTBDvsWT_scruff_volcano_data$diffexpressed[MTBDvsWT_scruff_volcano_data$L2FC > 0 & 
                                               MTBDvsWT_scruff_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < 0 and -log10pval > 1.3, set as "DOWN"
MTBDvsWT_scruff_volcano_data$diffexpressed[MTBDvsWT_scruff_volcano_data$L2FC < 0 & 
                                               MTBDvsWT_scruff_volcano_data$log10pval > 1.3] <- "DOWN"

# Add column used to label points
MTBDvsWT_scruff_volcano_data$DCMlabel <- ifelse(grepl("UP", MTBDvsWT_scruff_volcano_data$diffexpressed), 
                                                  MTBDvsWT_scruff_volcano_data$Metabolite, 
                                                  ifelse(grepl("DOWN", MTBDvsWT_scruff_volcano_data$diffexpressed), 
                                                         MTBDvsWT_scruff_volcano_data$Metabolite, NA))

volcano_fun(MTBDvsWT_scruff_volcano_data)
volcano_fun_label(MTBDvsWT_scruff_volcano_data)

## volcano_data plot for mdx post-scruff vs baseline

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
mdx_scruffvsbaseline_volcano_data$diffexpressed <- "NO"
# if log2Foldchange > 0 and -log10pval > 1.3, set as "UP"
mdx_scruffvsbaseline_volcano_data$diffexpressed[mdx_scruffvsbaseline_volcano_data$L2FC > 0 & 
                                                  mdx_scruffvsbaseline_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < 0 and -log10pval > 1.3, set as "DOWN"
mdx_scruffvsbaseline_volcano_data$diffexpressed[mdx_scruffvsbaseline_volcano_data$L2FC < 0 & 
                                             mdx_scruffvsbaseline_volcano_data$log10pval > 1.3] <- "DOWN"

# Add column used to label points
mdx_scruffvsbaseline_volcano_data$DCMlabel <- ifelse(grepl("UP", mdx_scruffvsbaseline_volcano_data$diffexpressed), 
                                                     mdx_scruffvsbaseline_volcano_data$Metabolite, 
                                                ifelse(grepl("DOWN", mdx_scruffvsbaseline_volcano_data$diffexpressed), 
                                                       mdx_scruffvsbaseline_volcano_data$Metabolite, NA))

volcano_fun(mdx_scruffvsbaseline_volcano_data)
volcano_fun_label(mdx_scruffvsbaseline_volcano_data)

## volcano_data plot for WT post-scruff vs baseline

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
WT_scruffvsbaseline_volcano_data$diffexpressed <- "NO"
# if log2Foldchange > 0 and -log10pval > 1.3, set as "UP"
WT_scruffvsbaseline_volcano_data$diffexpressed[WT_scruffvsbaseline_volcano_data$L2FC > 0 & 
                                                  WT_scruffvsbaseline_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < 0 and -log10pval > 1.3, set as "DOWN"
WT_scruffvsbaseline_volcano_data$diffexpressed[WT_scruffvsbaseline_volcano_data$L2FC < 0 & 
                                             WT_scruffvsbaseline_volcano_data$log10pval > 1.3] <- "DOWN"

# Add column used to label points
WT_scruffvsbaseline_volcano_data$DCMlabel <- ifelse(grepl("UP", WT_scruffvsbaseline_volcano_data$diffexpressed), 
                                                     WT_scruffvsbaseline_volcano_data$Metabolite, 
                                                     ifelse(grepl("DOWN", WT_scruffvsbaseline_volcano_data$diffexpressed), 
                                                            WT_scruffvsbaseline_volcano_data$Metabolite, NA))

volcano_fun(WT_scruffvsbaseline_volcano_data)
volcano_fun_label(WT_scruffvsbaseline_volcano_data)

## volcano_data plot for MTBD post-scruff vs baseline

# Add a column to the data frame to specify if they are UP- or DOWN- regulated
MTBD_scruffvsbaseline_volcano_data$diffexpressed <- "NO"
# if log2Foldchange > 0 and -log10pval > 1.3, set as "UP"
MTBD_scruffvsbaseline_volcano_data$diffexpressed[MTBD_scruffvsbaseline_volcano_data$L2FC > 0 & 
                                                 MTBD_scruffvsbaseline_volcano_data$log10pval > 1.3] <- "UP"
# if log2Foldchange < 0 and -log10pval > 1.3, set as "DOWN"
MTBD_scruffvsbaseline_volcano_data$diffexpressed[MTBD_scruffvsbaseline_volcano_data$L2FC < 0 & 
                                                 MTBD_scruffvsbaseline_volcano_data$log10pval > 1.3] <- "DOWN"

# Add column used to label points
MTBD_scruffvsbaseline_volcano_data$DCMlabel <- ifelse(grepl("UP", MTBD_scruffvsbaseline_volcano_data$diffexpressed), 
                                                    MTBD_scruffvsbaseline_volcano_data$Metabolite, 
                                                    ifelse(grepl("DOWN", MTBD_scruffvsbaseline_volcano_data$diffexpressed), 
                                                           MTBD_scruffvsbaseline_volcano_data$Metabolite, NA))

volcano_fun(MTBD_scruffvsbaseline_volcano_data)
volcano_fun_label(MTBD_scruffvsbaseline_volcano_data)

-----

# PLS-DA #####

# PLS discriminant analysis with automatic selection of components
data("iris")
my_pls = plsDA(iris[,1:4], iris$Species, autosel=TRUE)
my_pls$confusion
my_pls$error_rate

