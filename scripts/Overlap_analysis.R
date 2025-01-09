library(dplyr)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(Giotto)
library(reticulate)
library(circlize)
library(data.table)

use_python('/home/carlo/.local/share/r-miniconda/envs/giotto_env/bin/python', required = TRUE)

# Load data files
load("~/TumorGrowth/Wanglong/GSMobjfinal.Robj")
load("~/TumorGrowth/Wanglong/signfinal.Robj")
sign$EMR <- NULL

# Load Valdolivas datasets
V573_1 <- Load10X_Spatial('~/TumorGrowth/Valdolivas/SN048_A121573_Rep1')
V573_2 <- Load10X_Spatial('~/TumorGrowth/Valdolivas/SN048_A121573_Rep2')
V371_1 <- Load10X_Spatial('~/TumorGrowth/Valdolivas/SN048_A416371_Rep1')
V371_2 <- Load10X_Spatial('~/TumorGrowth/Valdolivas/SN048_A416371_Rep2')
V838_1 <- Load10X_Spatial('~/TumorGrowth/Valdolivas/SN84_A120838_Rep1')
V838_2 <- Load10X_Spatial('~/TumorGrowth/Valdolivas/SN84_A120838_Rep2')
V763_1 <- Load10X_Spatial('~/TumorGrowth/Valdolivas/SN123_A551763_Rep1')
V763_2 <- Load10X_Spatial('~/TumorGrowth/Valdolivas/SN124_A551763_Rep2')
V688_1 <- Load10X_Spatial('~/TumorGrowth/Valdolivas/SN123_A595688_Rep1')
V688_2 <- Load10X_Spatial('~/TumorGrowth/Valdolivas/SN124_A595688_Rep2')
V015_1 <- Load10X_Spatial('~/TumorGrowth/Valdolivas/SN123_A798015_Rep1')
V015_2 <- Load10X_Spatial('~/TumorGrowth/Valdolivas/SN124_A798015_Rep2')
V797_1 <- Load10X_Spatial('~/TumorGrowth/Valdolivas/SN123_A938797_Rep1_X')
V797_2 <- Load10X_Spatial('~/TumorGrowth/Valdolivas/SN124_A938797_Rep2')

# Custom subsetting function for Visium objects
custom_subset_visium <- function(seurat_object, cells) {
  # Subset the counts matrix
  subset_counts <- seurat_object@assays$Spatial@layers$counts[, cells, drop = FALSE]
  # Create a new Seurat object with the subsetted counts matrix and metadata
  new_seurat_object <- CreateSeuratObject(counts = subset_counts, meta.data = seurat_object@meta.data[cells, ])

  # Subset the image data if present
  if (!is.null(seurat_object@images)) {
    new_seurat_object@images <- seurat_object@images
    for (image_name in names(seurat_object@images)) {
      new_seurat_object@images[[image_name]]@boundaries$centroids@coords <- seurat_object@images[[image_name]]@boundaries$centroids@coords[cells, , drop = FALSE]
    }
  }

  # Subset the graph data if present
  if (!is.null(seurat_object@graphs)) {
    graph_object <- as.Graph(subset_counts)
    new_seurat_object@graphs <- list(graph = graph_object)
  }

  return(new_seurat_object)
}

# Create combined spatial object
allspatial <- CreateSeuratObject(counts = Matrix(cbind(GSM7058756_C1@assays$Spatial@layers$counts,
                                                     GSM7058757_C2@assays$Spatial@layers$counts,
                                                     GSM7058758_C3@assays$Spatial@layers$counts,
                                                     GSM7058759_C4@assays$Spatial@layers$counts), sparse = TRUE),
                               meta.data = as.data.frame(rbind(GSM7058756_C1@meta.data,
                                                               GSM7058757_C2@meta.data,
                                                               GSM7058758_C3@meta.data,
                                                               GSM7058759_C4@meta.data)))

spotstoconsider <- intersect(rownames(allspatial@meta.data)[allspatial$Fibroblasts_C2L > quantile(allspatial$Fibroblasts_C2L, .85)],
                             rownames(allspatial@meta.data)[allspatial$Fibroblasts_C2L_percentage > 20])

GSM7058759_C4_fibroenriched <- custom_subset_visium(GSM7058759_C4, intersect(rownames(GSM7058759_C4@meta.data), spotstoconsider))
GSM7058756_C1_fibroenriched <- custom_subset_visium(GSM7058756_C1, intersect(rownames(GSM7058756_C1@meta.data), spotstoconsider))
GSM7058757_C2_fibroenriched <- custom_subset_visium(GSM7058757_C2, intersect(rownames(GSM7058757_C2@meta.data), spotstoconsider))
GSM7058758_C3_fibroenriched <- custom_subset_visium(GSM7058758_C3, intersect(rownames(GSM7058758_C3@meta.data), spotstoconsider))
rm(allspatial)

# Function to process metadata for specific assay
deconv_metadata <- function(data, deconv_data, exp_id, cols, assay = "Spatial") {
  md <- deconv_data %>% filter(ST.exp == exp_id)
  md <- md[, cols]
  md <- md[md$ST.barcode %in% rownames(data@meta.data), ]
  data <- data[, md$ST.barcode]
  rownames(md) <- rownames(data@meta.data)
  colnames(md) <- gsub('$', '_Deconv', colnames(md))
  data@meta.data <- cbind(data@meta.data, md)
  return(data)
}

# Function to process C2L for Korean deconvolution
process_C2L_Korean <- function(seurat_object, df = FALSE) {
  if (df == FALSE) {
    mt <- seurat_object@meta.data
    cell_type_columns <- mt[, stringr::str_detect(pattern = "^(?!ST\\.exp_)(?!ST\\.barcode_).+_KoreanDeconv$", string = colnames(mt))]
    row_totals <- rowSums(cell_type_columns)
    cell_type_percentages <- sweep(cell_type_columns, 1, row_totals, "/") * 100
    colnames(cell_type_percentages) <- paste0(colnames(cell_type_percentages), "_percentageKorean")
    mt <- cbind(mt, cell_type_percentages)
    seurat_object@meta.data <- mt
    return(seurat_object)
  } else {
    cell_type_columns <- seurat_object[, stringr::str_detect(pattern = "^(?!spot.id_)(?!ST\\.exp_)(?!ST\\.barcode_).+_KoreanDeconv$", string = colnames(seurat_object))]
    row_totals <- rowSums(cell_type_columns)
    cell_type_percentages <- sweep(cell_type_columns, 1, row_totals, "/") * 100
    colnames(cell_type_percentages) <- paste0(colnames(cell_type_percentages), "_percentageKorean")
    seurat_object <- cbind(seurat_object, cell_type_percentages)
    return(seurat_object)
  }
}

# Function for PAGE analysis
process_PAGEanalysis <- function(seurat_object, signatures_all, only_fibro = TRUE, assay = "Spatial") {
  if (only_fibro) {
    signatures_all <- signatures_all[!names(signatures_all) %in% c("CCM", "EMRM")]
  }

  raw_exprs <- if (assay == "Spatial") {
    seurat_object@assays$Spatial@layers$counts
  } else if (assay == "RNA") {
    seurat_object@assays$RNA@counts
  } else {
    stop("Unsupported assay type")
  }

  rownames(raw_exprs) <- rownames(seurat_object@assays[[assay]]@features@.Data)
  colnames(raw_exprs) <- rownames(seurat_object@meta.data)
  spatial_locs <- as.data.table(GetTissueCoordinates(seurat_object)[, 1:2])
  colnames(spatial_locs) <- c("x", "y")

  myGiottoObj <- createGiottoObject(raw_exprs = raw_exprs, spatial_locs = spatial_locs)
  myGiottoObj <- normalizeGiotto(gobject = myGiottoObj)

  all_signatures <- names(signatures_all)
  signature_matrix_complete <- makeSignMatrixPAGE(
    sign_names = all_signatures,
    sign_list = lapply(all_signatures, function(sign) signatures_all[[sign]])
  )

  myGiottoObj_initial <- runPAGEEnrich(
    gobject = myGiottoObj,
    sign_matrix = signature_matrix_complete,
    min_overlap_genes = 2,
    output_enrichment = c("original", "zscore")
  )

  myGiottoObj_pval <- runPAGEEnrich(
    gobject = myGiottoObj,
    sign_matrix = signature_matrix_complete,
    min_overlap_genes = 2,
    p_value = TRUE,
    reverse_log_scale = FALSE
  )

  pval_df <- as.data.frame(myGiottoObj_pval@spatial_enrichment$PAGE)
  zscore_df <- as.data.frame(myGiottoObj_initial@spatial_enrichment$PAGE)
  colnames(zscore_df) <- paste0(colnames(zscore_df), "PAGE")
  seurat_object@meta.data <- cbind(seurat_object@meta.data, zscore_df)

  significant_spots <- pval_df > 1.301
  significant_df <- as.data.frame(significant_spots)
  colnames(significant_df) <- paste0(colnames(significant_df), "PAGE_significant")
  seurat_object@meta.data <- cbind(seurat_object@meta.data, significant_df)

  return(seurat_object)
}

# Function to extract metadata by columns
extract_metadata <- function(seurat_object, start_col, end_col) {
  return(seurat_object@meta.data[, start_col:end_col])
}

# Fisher's z-transformation and its inverse
fisher_z <- function(r) {
  return(0.5 * log((1 + r) / (1 - r)))
}

inverse_fisher_z <- function(z) {
  return((exp(2 * z) - 1) / (exp(2 * z) + 1))
}

# Example: Extract metadata and calculate correlations
metadata_list <- list(
  extract_metadata(V763_1_fibro, 139, 156),
  extract_metadata(V371_2_fibro, 181, 198),
  extract_metadata(V688_1_fibro, 139, 156),
  extract_metadata(V015_1_fibro, 139, 156),
  extract_metadata(V015_2_fibro, 139, 156),
  extract_metadata(V797_1_fibro, 139, 156),
  extract_metadata(V797_2_fibro, 139, 156),
  extract_metadata(V838_1_fibro, 139, 156),
  extract_metadata(V838_2_fibro, 139, 156),
  extract_metadata(GSM7058756_C1_fibroenriched, 77, 94),
  extract_metadata(GSM7058757_C2_fibroenriched, 77, 94),
  extract_metadata(GSM7058759_C4_fibroenriched, 77, 94)
)

correlation_matrices <- lapply(metadata_list, cor)
num_spots <- sapply(metadata_list, nrow)

# Fisher's z-transformation to each matrix
z_matrices <- lapply(correlation_matrices, function(mat) {
  apply(mat, c(1, 2), fisher_z)
})

# Weight the transformed values by the number of spots
weighted_z_matrices <- mapply(function(z_matrix, n_spots) {
  z_matrix * n_spots
}, z_matrices, num_spots, SIMPLIFY = FALSE)

# Compute the sum of weighted z-values and total number of spots
sum_weighted_z_matrix <- Reduce("+", weighted_z_matrices)
total_spots <- sum(num_spots)

# Compute the weighted average of z-values
mean_weighted_z_matrix <- sum_weighted_z_matrix / total_spots

# Apply the inverse Fisher's z-transformation
merged_correlation_matrix <- apply(mean_weighted_z_matrix, c(1, 2), inverse_fisher_z)

# Ensure the diagonal remains 1 (since it's a correlation matrix)
diag(merged_correlation_matrix) <- 1

# Perform hierarchical clustering
row_model <- hclust(dist(merged_correlation_matrix))
col_model <- hclust(t(dist(merged_correlation_matrix)))

library(ComplexHeatmap)
# Create a Heatmap object
ht <- Heatmap(merged_correlation_matrix, 
              name = "Merged\nCorrelation\ncorrected for n.spots",
              col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
              cluster_rows = row_model,
              cluster_columns = col_model,
              show_row_dend = TRUE,
              show_column_dend = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_dend_side = "left",
              column_dend_side = "top",
              row_names_side = "left",
              column_names_side = "bottom",
              heatmap_legend_param = list(title = "Merged\nCorrelation\ncorrected\nfor n.spots"),
              row_dend_width = unit(4, "cm"),
              column_dend_height = unit(4, "cm"))

# Save the heatmap
Cairo::CairoPNG("FinalMergedCorr_matrix_85_fixed20_final.png", width = 1000, height = 1000)
draw(ht, heatmap_legend_side = "right")
dev.off()

# Function to calculate hypergeometric p-values
calculate_hypergeometric_pvals <- function(df) {
  categories <- colnames(df)
  comb_matrix <- combn(categories, 2, simplify = FALSE)

  pval_matrix <- matrix(0, nrow = length(categories), ncol = length(categories))
  colnames(pval_matrix) <- categories
  rownames(pval_matrix) <- categories

  N <- nrow(df)

  for (comb in comb_matrix) {
    cat1 <- comb[1]
    cat2 <- comb[2]

    K <- sum(df[[cat1]])
    n <- sum(df[[cat2]])
    k <- sum(df[[cat1]] & df[[cat2]])

    pval <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)

    pval_matrix[cat1, cat2] <- pval
    pval_matrix[cat2, cat1] <- pval
  }

  diag(pval_matrix) <- 1

  return(pval_matrix)
}

# Example: Calculate hypergeometric p-values for fibro datasets
pval_matrices <- lapply(metadata_list, calculate_hypergeometric_pvals)

# Function to combine p-values using Fisher's method
combine_pvals <- function(pval_matrices) {
  categories <- rownames(pval_matrices[[1]])
  combined_pval_matrix <- matrix(0, nrow = length(categories), ncol = length(categories))
  colnames(combined_pval_matrix) <- categories
  rownames(combined_pval_matrix) <- categories

  for (i in 1:length(categories)) {
    for (j in 1:length(categories)) {
      pvals <- sapply(pval_matrices, function(x) x[i, j])
      if (all(pvals == 1)) {
        combined_pval <- 1
      } else {
        chisq_stat <- -2 * sum(log(pvals))
        combined_pval <- pchisq(chisq_stat, df = 2 * length(pvals), lower.tail = FALSE)
      }
      combined_pval_matrix[i, j] <- combined_pval
    }
  }

  return(combined_pval_matrix)
}

# Combine p-values from all datasets
combined_pval_matrix <- combine_pvals(pval_matrices)

# Perform hierarchical clustering on p-values
row_model <- hclust(dist(combined_pval_matrix))
col_model <- hclust(t(dist(combined_pval_matrix)))

# Create a heatmap for combined p-values - we are adding a minimal value to the matrix to avoid log(0)
ht <- Heatmap(-log10(combined_pval_matrix + 5e-324), 
              name = "Combined\nenrichment\nscore\nvia Fisher's\nMethod",
              col = colorRamp2(c(0, 20, 400), c("white", "#E94B3C", "#6C1413")),
              cluster_rows = row_model,
              cluster_columns = col_model,
              show_row_dend = TRUE,
              show_column_dend = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              row_dend_side = "left",
              column_dend_side = "top",
              row_names_side = "left",
              column_names_side = "bottom",
              heatmap_legend_param = list(title = "Combined\nenrichment\nscore\nvia Fisher's\nMethod"),
              row_dend_width = unit(3, "cm"),
              column_dend_height = unit(4, "cm"))

# Save the heatmap
Cairo::CairoPNG("Hypergeometri_p_values_F_cmplt_nomrCAF_completelink.png", width = 1000, height = 1000)
draw(ht, heatmap_legend_side = "right")
dev.off()
