

######################################################################################## Percentile Heatmap between Healthy, Uninflamed and Inflamed (figure 6B)

library(qs)


emr_new<-read.csv('EMR_New.csv')
mrCAF_new<-read.csv('mrCAF_New.csv')

setwd('/media/carlo/FFF9-F0A1/NatureCancer_2024/Fibroblasts/')

fibro_gse<-mrCAF_new%>%filter(Dataset%in%c('C.Smillie','G.Li','GSE125527'))

fibro_gse<-mrCAF_new

fibro_gse <- fibro_gse %>%
  mutate(Patient_Class = paste(Class, Patient, sep = "_"))

fibro_gse<-fibro_gse%>%filter(Class%in%c('Healthy','Uninflamed','Inflamed'))
rm(mrCAF_new)


emr_malign<-emr_new
emr_malign <- emr_malign %>%
  mutate(Patient_Class = paste(Class, Patient, sep = "_"))
emr_malign<-emr_new%>%filter(Dataset%in%c('C.Smillie','G.Li','GSE125527'))
emr_malign<-emr_malign%>%filter(Class%in%c('Healthy','Uninflamed','Inflamed'))
rm(emr_new)

setwd('/media/carlo/FFF9-F0A1/NatureCancer_2024/Fibroblasts/')
fibro_mtx<-qread('SparseMatrix_fibro.qs')

setwd('/media/carlo/FFF9-F0A1/NatureCancer_2024/')
epi_mtx <- qread('AllMalign_sparse.qs')

epi_mtx <- epi_mtx[,emr_malign$X.1]

HealthyInflamedUninflamed_mtx <- fibro_mtx[fibro_gse$X.1,]


HealthyInflamedUninflamed <- CreateSeuratObject(counts = t(HealthyInflamedUninflamed_mtx), meta.data = fibro_gse)
Idents(HealthyInflamedUninflamed) = 'Class'

scores <- HealthyInflamedUninflamed@meta.data$mrCAFMean

# Define quantile thresholds from 75% to 95% by 1% increments and add 99%
quantiles_seq <- seq(0.75, 0.95, by = 0.01)
quantiles_seq <- c(quantiles_seq, 0.99)
thresh_values <- quantile(scores, quantiles_seq)
print(thresh_values)

# Initialize lists to store output
prop_list <- list()
odds_list_Healthy <- list()
odds_list_Uninflamed <- list()
odds_list_NotInflamed <- list()
pval_list_Healthy <- list()
pval_list_Uninflamed <- list()
pval_list_NotInflamed <- list()

# Loop over each threshold
for(i in seq_along(thresh_values)){
  threshold_value <- thresh_values[i]
  thresh_label <- paste0(round(quantiles_seq[i] * 100), "%")
  
  # Define a column name for the binary positivity
  pos_colname <- paste0("mrCAF_pos_", thresh_label)
  
  # Add a binary metadata column: TRUE if the cell's score > threshold
  HealthyInflamedUninflamed <- AddMetaData(
    HealthyInflamedUninflamed,
    metadata = scores > threshold_value,
    col.name = pos_colname
  )
  
  # Extract updated metadata
  meta <- HealthyInflamedUninflamed@meta.data
  
  # Build a contingency table: rows = Class, columns = positivity (TRUE/FALSE)
  tab <- table(meta$Class, meta[[pos_colname]])
  
  # Compute proportions for each class
  prop_H <- tab["Healthy", "TRUE"] / sum(tab["Healthy", ])
  prop_I <- tab["Inflamed", "TRUE"] / sum(tab["Inflamed", ])
  prop_U <- tab["Uninflamed", "TRUE"] / sum(tab["Uninflamed", ])
  prop_list[[thresh_label]] <- c(Healthy = prop_H, Inflamed = prop_I, Uninflamed = prop_U)
  
  # --- Manual Odds Ratio Calculation ---
  # Extract counts for each class:
  # For Inflamed:
  I_true <- tab["Inflamed", "TRUE"]
  I_false <- tab["Inflamed", "FALSE"]
  # For Healthy:
  H_true <- tab["Healthy", "TRUE"]
  H_false <- tab["Healthy", "FALSE"]
  # For Uninflamed:
  U_true <- tab["Uninflamed", "TRUE"]
  U_false <- tab["Uninflamed", "FALSE"]
  
  # Odds Ratio: Inflamed vs Healthy
  if(I_false == 0 || H_false == 0){
    OR_IH <- NA
  } else {
    OR_IH <- (I_true / I_false) / (H_true / H_false)
  }
  odds_list_Healthy[[thresh_label]] <- OR_IH
  # p-value for Inflamed vs Healthy using fisher.test:
  tab_IH <- rbind(Inflamed = tab["Inflamed", ], Healthy = tab["Healthy", ])
  fisher_IH <- fisher.test(tab_IH)
  pval_list_Healthy[[thresh_label]] <- fisher_IH$p.value
  
  # Odds Ratio: Inflamed vs Uninflamed
  if(I_false == 0 || U_false == 0){
    OR_IU <- NA
  } else {
    OR_IU <- (I_true / I_false) / (U_true / U_false)
  }
  odds_list_Uninflamed[[thresh_label]] <- OR_IU
  # p-value for Inflamed vs Uninflamed:
  tab_IU <- rbind(Inflamed = tab["Inflamed", ], Uninflamed = tab["Uninflamed", ])
  fisher_IU <- fisher.test(tab_IU)
  pval_list_Uninflamed[[thresh_label]] <- fisher_IU$p.value
  
  # Odds Ratio: Inflamed vs NotInflamed (union of Healthy and Uninflamed)
  NotInflamed_true <- H_true + U_true
  NotInflamed_false <- H_false + U_false
  if(I_false == 0 || NotInflamed_false == 0){
    OR_INI <- NA
  } else {
    OR_INI <- (I_true / I_false) / (NotInflamed_true / NotInflamed_false)
  }
  odds_list_NotInflamed[[thresh_label]] <- OR_INI
  # p-value for Inflamed vs NotInflamed:
  notinflamed_counts <- colSums(tab[c("Healthy", "Uninflamed"), ])
  tab_INI <- rbind(Inflamed = tab["Inflamed", ], NotInflamed = notinflated_counts)
  # Correction: use the computed notinflamed_counts
  tab_INI <- rbind(Inflamed = tab["Inflamed", ], NotInflamed = notinflamed_counts)
  fisher_INI <- fisher.test(tab_INI)
  pval_list_NotInflamed[[thresh_label]] <- fisher_INI$p.value
}

# Combine proportions into a data frame for visualization (if needed)
heatmap_df <- do.call(rbind, prop_list)
heatmap_df <- as.data.frame(heatmap_df)
heatmap_df$Threshold <- rownames(heatmap_df)
df_tidy <- melt(heatmap_df, id.vars = "Threshold",
                variable.name = "Class", value.name = "Proportion")
df_tidy$Threshold <- factor(df_tidy$Threshold, levels = names(thresh_values))
df_tidy$Class <- factor(df_tidy$Class, levels = c("Healthy", "Uninflamed", "Inflamed"))

# Create data frames for the manual odds ratios
stats_df_Healthy <- data.frame(
  Threshold = names(odds_list_Healthy),
  p_value = unlist(pval_list_Healthy),
  odds_ratio = unlist(odds_list_Healthy),
  Comparison = "Inflamed vs Healthy"
)
stats_df_Uninflamed <- data.frame(
  Threshold = names(odds_list_Uninflamed),
  p_value = unlist(pval_list_Uninflamed),
  odds_ratio = unlist(odds_list_Uninflamed),
  Comparison = "Inflamed vs Uninflamed"
)
stats_df_NotInflamed <- data.frame(
  Threshold = names(odds_list_NotInflamed),
  p_value = unlist(pval_list_NotInflamed),
  odds_ratio = unlist(odds_list_NotInflamed),
  Comparison = "Inflamed vs NotInflamed"
)

# Combine all statistics into one table
stats_all <- rbind(stats_df_Healthy, stats_df_Uninflamed, stats_df_NotInflamed)
stats_all$Threshold <- factor(stats_all$Threshold, levels = names(thresh_values))

write.csv(stats_all, "Supplementary_Stats_mrCAF.csv", row.names = FALSE)

# Create the heatmap using ggplot2
p <- ggplot(df_tidy, aes(x = Class, y = Threshold, fill = Proportion)) +
  geom_tile(color = "white") +
  # Using scale_fill_gradient2 to set a custom midpoint (0.75 here)
  scale_fill_gradient2(low = "white", high = "darkred", 
                       midpoint = 0.1, name = "Fraction of \nmrCAF+ cells") +
  theme_minimal() + 
  labs(x = NULL, y = "Thresholds", title = NULL, subtitle = NULL)+
  # Remove titles and axis labels for a cleaner look
  labs(x = NULL, y = NULL, title = NULL, subtitle = NULL) +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 12),
        panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot to PDF
pdf('ggpubr_Inflamed_vs_Healthy_Uninflamed_5thQuantiles.pdf', width = 3.5, height = 10)
print(p)
dev.off()


Idents(epi_mtx) = 'Dataset'

epi_mtx <- subset(epi_mtx,idents = c('C.Smillie',
                                     'G.Li',
                                     'GSE125527'))

scores <- epi_mtx@meta.data$EMRMean

# Define quantile thresholds from 75% to 95% by 1% increments and add 99%
quantiles_seq <- seq(0.75, 0.95, by = 0.01)
quantiles_seq <- c(quantiles_seq, 0.99)
thresh_values <- quantile(scores, quantiles_seq)
print(thresh_values)

# Initialize lists to store output
prop_list <- list()
odds_list_Healthy <- list()
odds_list_Uninflamed <- list()
odds_list_NotInflamed <- list()
pval_list_Healthy <- list()
pval_list_Uninflamed <- list()
pval_list_NotInflamed <- list()

# Loop over each threshold
for(i in seq_along(thresh_values)){
  threshold_value <- thresh_values[i]
  thresh_label <- paste0(round(quantiles_seq[i] * 100), "%")
  
  # Define a column name for the binary positivity
  pos_colname <- paste0("EMR_pos_", thresh_label)
  
  # Add a binary metadata column: TRUE if the cell's score > threshold
  HealthyInflamedUninflamed <- AddMetaData(
    HealthyInflamedUninflamed,
    metadata = scores > threshold_value,
    col.name = pos_colname
  )
  
  # Extract updated metadata
  meta <- HealthyInflamedUninflamed@meta.data
  
  # Build a contingency table: rows = Class, columns = positivity (TRUE/FALSE)
  tab <- table(meta$Class, meta[[pos_colname]])
  
  # Compute proportions for each class
  prop_H <- tab["Healthy", "TRUE"] / sum(tab["Healthy", ])
  prop_I <- tab["Inflamed", "TRUE"] / sum(tab["Inflamed", ])
  prop_U <- tab["Uninflamed", "TRUE"] / sum(tab["Uninflamed", ])
  prop_list[[thresh_label]] <- c(Healthy = prop_H, Inflamed = prop_I, Uninflamed = prop_U)
  
  # --- Manual Odds Ratio Calculation ---
  # Extract counts for each class:
  # For Inflamed:
  I_true <- tab["Inflamed", "TRUE"]
  I_false <- tab["Inflamed", "FALSE"]
  # For Healthy:
  H_true <- tab["Healthy", "TRUE"]
  H_false <- tab["Healthy", "FALSE"]
  # For Uninflamed:
  U_true <- tab["Uninflamed", "TRUE"]
  U_false <- tab["Uninflamed", "FALSE"]
  
  # Odds Ratio: Inflamed vs Healthy
  if(I_false == 0 || H_false == 0){
    OR_IH <- NA
  } else {
    OR_IH <- (I_true / I_false) / (H_true / H_false)
  }
  odds_list_Healthy[[thresh_label]] <- OR_IH
  # p-value for Inflamed vs Healthy using fisher.test:
  tab_IH <- rbind(Inflamed = tab["Inflamed", ], Healthy = tab["Healthy", ])
  fisher_IH <- fisher.test(tab_IH)
  pval_list_Healthy[[thresh_label]] <- fisher_IH$p.value
  
  # Odds Ratio: Inflamed vs Uninflamed
  if(I_false == 0 || U_false == 0){
    OR_IU <- NA
  } else {
    OR_IU <- (I_true / I_false) / (U_true / U_false)
  }
  odds_list_Uninflamed[[thresh_label]] <- OR_IU
  # p-value for Inflamed vs Uninflamed:
  tab_IU <- rbind(Inflamed = tab["Inflamed", ], Uninflamed = tab["Uninflamed", ])
  fisher_IU <- fisher.test(tab_IU)
  pval_list_Uninflamed[[thresh_label]] <- fisher_IU$p.value
  
  # Odds Ratio: Inflamed vs NotInflamed (union of Healthy and Uninflamed)
  NotInflamed_true <- H_true + U_true
  NotInflamed_false <- H_false + U_false
  if(I_false == 0 || NotInflamed_false == 0){
    OR_INI <- NA
  } else {
    OR_INI <- (I_true / I_false) / (NotInflamed_true / NotInflamed_false)
  }
  odds_list_NotInflamed[[thresh_label]] <- OR_INI
  # p-value for Inflamed vs NotInflamed:
  notinflamed_counts <- colSums(tab[c("Healthy", "Uninflamed"), ])
  tab_INI <- rbind(Inflamed = tab["Inflamed", ], NotInflamed = notinflamed_counts)
  # Correction: use the computed notinflamed_counts
  tab_INI <- rbind(Inflamed = tab["Inflamed", ], NotInflamed = notinflamed_counts)
  fisher_INI <- fisher.test(tab_INI)
  pval_list_NotInflamed[[thresh_label]] <- fisher_INI$p.value
}

# Combine proportions into a data frame for visualization (if needed)
heatmap_df <- do.call(rbind, prop_list)
heatmap_df <- as.data.frame(heatmap_df)
heatmap_df$Threshold <- rownames(heatmap_df)
df_tidy <- melt(heatmap_df, id.vars = "Threshold",
                variable.name = "Class", value.name = "Proportion")
df_tidy$Threshold <- factor(df_tidy$Threshold, levels = names(thresh_values))
df_tidy$Class <- factor(df_tidy$Class, levels = c("Healthy", "Uninflamed", "Inflamed"))

# Create data frames for the manual odds ratios
stats_df_Healthy <- data.frame(
  Threshold = names(odds_list_Healthy),
  p_value = unlist(pval_list_Healthy),
  odds_ratio = unlist(odds_list_Healthy),
  Comparison = "Inflamed vs Healthy"
)
stats_df_Uninflamed <- data.frame(
  Threshold = names(odds_list_Uninflamed),
  p_value = unlist(pval_list_Uninflamed),
  odds_ratio = unlist(odds_list_Uninflamed),
  Comparison = "Inflamed vs Uninflamed"
)
stats_df_NotInflamed <- data.frame(
  Threshold = names(odds_list_NotInflamed),
  p_value = unlist(pval_list_NotInflamed),
  odds_ratio = unlist(odds_list_NotInflamed),
  Comparison = "Inflamed vs NotInflamed"
)

# Combine all statistics into one table
stats_all <- rbind(stats_df_Healthy, stats_df_Uninflamed, stats_df_NotInflamed)
stats_all$Threshold <- factor(stats_all$Threshold, levels = names(thresh_values))

write.csv(stats_all, "Supplementary_Stats_EMR.csv", row.names = FALSE)

# Create the heatmap using ggplot2
p <- ggplot(df_tidy, aes(x = Class, y = Threshold, fill = Proportion)) +
  geom_tile(color = "white") +
  # Using scale_fill_gradient2 to set a custom midpoint (0.75 here)
  scale_fill_gradient2(low = "white", high = "darkgreen", 
                       midpoint = 0.1, name = "Fraction of \nEMR+ cells") +
  theme_minimal() + 
  labs(x = NULL, y = "Thresholds", title = NULL, subtitle = NULL)+
  # Remove titles and axis labels for a cleaner look
  labs(x = NULL, y = NULL, title = NULL, subtitle = NULL) +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 12),
        panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot to PDF
pdf('ggpubr_Inflamed_vs_Healthy_Uninflamed_5thQuantiles_EMR.pdf', width = 3.5, height = 10)
print(p)
dev.off()
