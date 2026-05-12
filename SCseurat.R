# Single Cell RNA-seq Analysis using Seurat
# Sample: GSM8615073_302

library(Seurat)
library(tidyverse)

# Settings --------------------------------------------------------------------

project_dir <- "/Users/ghazl/Projects/Charite_Internship"
sample_id <- "GSM8615073_302"

data_dir <- file.path(project_dir, "GSE281219_10x", sample_id)
output_dir <- file.path(project_dir, "R")

min_cells <- 3
min_features_create <- 200
min_features_qc <- 200
max_features_qc <- 5000
max_mt_percent <- 20

variable_features <- 2000
pca_dimensions <- 30
analysis_dimensions <- 1:20
cluster_resolution <- 0.5

dir.create(output_dir, showWarnings = FALSE)


# Helper functions ------------------------------------------------------------

step_message <- function(step, message) {
  cat("\nSTEP ", step, ": ", message, "\n", sep = "")
}

done_message <- function(message) {
  cat("Done: ", message, "\n", sep = "")
}

output_file <- function(filename) {
  file.path(output_dir, filename)
}

save_png <- function(filename, plot, width = 800, height = 600) {
  png(output_file(filename), width = width, height = height)
  print(plot)
  dev.off()
}


# 1. Load data ----------------------------------------------------------------

step_message(1, "Loading data")

cat("Working directory: ", getwd(), "\n", sep = "")
cat("Data directory: ", data_dir, "\n", sep = "")
cat("Sample: ", sample_id, "\n", sep = "")

required_files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
missing_files <- required_files[
  !file.exists(file.path(data_dir, required_files))
]

if (length(missing_files) > 0) {
  stop(
    "Missing 10X files in ",
    data_dir,
    ": ",
    paste(missing_files, collapse = ", ")
  )
}

seurat_data <- Read10X(data.dir = data_dir)
total_cells <- ncol(seurat_data)
total_features <- nrow(seurat_data)

seurat_obj <- CreateSeuratObject(
  counts = seurat_data,
  project = sample_id,
  min.cells = min_cells,
  min.features = min_features_create
)

done_message("data loaded")
cat("Features: ", nrow(seurat_obj), "\n", sep = "")
cat("Cells: ", ncol(seurat_obj), "\n", sep = "")


# 2. Quality control ----------------------------------------------------------

step_message(2, "Calculating QC metrics")

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^MT-"
)

save_png(
  "01_qc_violin_plots.png",
  VlnPlot(
    seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3
  ),
  width = 1200,
  height = 400
)

done_message("QC metrics calculated")
cat(
  "Median features per cell: ",
  median(seurat_obj$nFeature_RNA),
  "\n",
  sep = ""
)
cat("Median counts per cell: ", median(seurat_obj$nCount_RNA), "\n", sep = "")
cat("Median percent MT: ", median(seurat_obj$percent.mt), "\n", sep = "")


# 3. Filter cells -------------------------------------------------------------

step_message(3, "Filtering low-quality cells")

seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > min_features_qc &
    nFeature_RNA < max_features_qc &
    percent.mt < max_mt_percent
)

done_message("filtering complete")
cat("Remaining cells: ", ncol(seurat_obj), "\n", sep = "")
cat("Remaining features: ", nrow(seurat_obj), "\n", sep = "")


# 4. Normalize and select features -------------------------------------------

step_message(4, "Normalizing data")

seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

done_message("data normalized")

step_message(5, "Finding highly variable genes")

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = variable_features
)

save_png(
  "05_variable_features_plot.png",
  VariableFeaturePlot(seurat_obj)
)

done_message("highly variable genes selected")


# 5. PCA, UMAP, and clustering ------------------------------------------------

step_message(6, "Scaling data")

seurat_obj <- ScaleData(seurat_obj)

done_message("data scaled")

step_message(7, "Running PCA")

seurat_obj <- RunPCA(
  seurat_obj,
  features = VariableFeatures(seurat_obj),
  npcs = pca_dimensions
)

save_png("04_pca_plot.png", DimPlot(seurat_obj, reduction = "pca"))
save_png(
  "06_pca_loadings.png",
  VizDimLoadings(seurat_obj, dims = 1:2, nfeatures = 15),
  width = 1000
)
save_png(
  "07_elbow_plot.png",
  ElbowPlot(seurat_obj, ndims = pca_dimensions)
)

done_message("PCA complete")

step_message(8, "Running UMAP")

seurat_obj <- RunUMAP(seurat_obj, dims = analysis_dimensions)

save_png("03_umap_plot.png", DimPlot(seurat_obj, reduction = "umap"))

done_message("UMAP complete")

step_message(9, "Finding neighbors and clusters")

seurat_obj <- FindNeighbors(seurat_obj, dims = analysis_dimensions)
seurat_obj <- FindClusters(seurat_obj, resolution = cluster_resolution)

save_png(
  "08_clusters_umap.png",
  DimPlot(seurat_obj, reduction = "umap", label = TRUE)
)

done_message("clustering complete")
cat("Number of clusters: ", length(unique(Idents(seurat_obj))), "\n", sep = "")


# 6. Marker genes -------------------------------------------------------------

step_message(10, "Finding cluster-specific markers")

markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  markers,
  output_file("cluster_markers.csv"),
  row.names = FALSE
)

done_message("markers identified")
cat("Total marker genes: ", nrow(markers), "\n", sep = "")
cat("\nTop 10 markers per cluster:\n")

markers |>
  group_by(cluster) |>
  slice_head(n = 10) |>
  print(n = Inf)

top_markers <- markers |>
  group_by(cluster) |>
  slice_head(n = 5) |>
  pull(gene)

save_png(
  "09_markers_heatmap.png",
  DoHeatmap(seurat_obj, features = unique(top_markers)),
  width = 1000
)

clusters_to_plot <- head(unique(markers$cluster), 3)

for (cluster_id in clusters_to_plot) {
  top_gene <- markers |>
    filter(cluster == cluster_id) |>
    slice_head(n = 1) |>
    pull(gene)

  save_png(
    paste0("10_marker_", cluster_id, "_", top_gene, ".png"),
    FeaturePlot(seurat_obj, features = top_gene)
  )
}


# 7. Save results -------------------------------------------------------------

step_message(11, "Saving analysis results")

saveRDS(
  seurat_obj,
  output_file(paste0("seurat_", sample_id, "_processed.rds"))
)

cluster_data <- data.frame(
  barcode = colnames(seurat_obj),
  cluster = Idents(seurat_obj)
)

write.csv(
  cluster_data,
  output_file(paste0("seurat_", sample_id, "_clusters.csv")),
  row.names = FALSE
)

qc_summary <- data.frame(
  metric = c(
    "Total cells",
    "Cells after QC",
    "Total features",
    "Mean features per cell",
    "Mean counts per cell",
    "Mean percent MT"
  ),
  value = c(
    format(total_cells, big.mark = ","),
    format(ncol(seurat_obj), big.mark = ","),
    format(total_features, big.mark = ","),
    round(mean(seurat_obj$nFeature_RNA), 2),
    round(mean(seurat_obj$nCount_RNA), 2),
    round(mean(seurat_obj$percent.mt), 2)
  )
)

write.table(
  qc_summary,
  output_file("qc_metrics_summary.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

done_message("analysis complete")
cat("Files saved to: ", output_dir, "\n", sep = "")
cat("\nMain output files:\n")
cat("- seurat_", sample_id, "_processed.rds\n", sep = "")
cat("- seurat_", sample_id, "_clusters.csv\n", sep = "")
cat("- cluster_markers.csv\n")
cat("- qc_metrics_summary.txt\n")
cat("- PNG visualization files\n")
