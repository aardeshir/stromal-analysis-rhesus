# Stromal Cell Analysis in Primate Colon
# Author: A. Ardeshir (Ardeshirlab - aardeshir@tulane.edu)
# Date: 2025-01-08

# Setup and Data Loading ----
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(SpatialExperiment)
  library(scater)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
  library(viridis)
  library(gridExtra)
})

# Create output directories
dir.create(file.path("output", "stromal_analysis", "plots"), 
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("output", "stromal_analysis", "processed"), 
           recursive = TRUE, showWarnings = FALSE)

# Load data
metadata_add <- read.csv(file.path("data", "Additional_metadata.csv"), header = TRUE)
gcm_PR03 <- readRDS(file.path("PR03_Chronic", "data", "processed", "counts.RDS"))
locs_PR03 <- readRDS(file.path("PR03_Chronic", "data", "processed", "xy.RDS"))
negcounts_PRO3 <- readRDS(file.path("PR03_Chronic", "data", "processed", "negcounts.RDS"))
metadata_PR03 <- readRDS(file.path("PR03_Chronic", "data", "processed", "metadata.RDS"))

# Process metadata and create SE object
metadata_PR03_merge <- merge(metadata_PR03, metadata_add)
metadata_PR03_merge$Run_Tissue_name <- NULL
rownames(metadata_PR03_merge) <- metadata_PR03_merge$cell_id
metadata_PR03_merge$negmean <- rowMeans(negcounts_PRO3)
cd_PR03 <- cbind(locs_PR03, metadata_PR03_merge)

se <- SpatialExperiment(
  assays = list(counts = t(gcm_PR03)),
  colData = cd_PR03
)

# Clean up
rm(gcm_PR03, locs_PR03, negcounts_PRO3, metadata_PR03, metadata_PR03_merge, cd_PR03, metadata_add)

# Normalize data
se <- computeLibraryFactors(se)
assay(se, "normcounts") <- normalizeCounts(se, log = FALSE)

# Subset to relevant animals and set proper condition labels
se <- se[, se$Animal_ID %in% c("RC56.1", "RB45"), drop = FALSE]
colData(se)$condition <- ifelse(colData(se)$Animal_ID == "RB45", "ICD", "Control")

# Define marker sets with technical justification
markers <- list(
  # Primary stromal markers - core identity
  stromal = c(
    "PDGFRA",  # Canonical stromal marker
    "NCAM1",   # Neural cell adhesion molecule
    "VIM"      # Mesenchymal marker
  ),
  # Negative selection - exclude other lineages
  negative = c(
    "PTPRC",   # CD45, exclude immune cells
    "PECAM1",  # Exclude endothelial cells
    "EPCAM"    # Exclude epithelial cells
  ),
  # Additional characterization
  additional = c(
    "CD9",     # Surface marker
    "CD34"     # Progenitor marker
  ),
  # Target gene
  target = "BMP2"
)

# Functions ----

# Function to analyze marker quality
analyze_marker_quality <- function(se, available_markers) {
  expr_mat <- assay(se, "normcounts")[available_markers,]
  
  # Calculate basic statistics
  marker_stats <- data.frame(
    Marker = available_markers,
    Detected_Cells = apply(expr_mat > 0, 1, sum),
    Mean_Expression = rowMeans(expr_mat),
    Median_Expression = apply(expr_mat, 1, median),
    Q75_Expression = apply(expr_mat, 1, function(x) quantile(x, 0.75)),
    Detection_Rate = apply(expr_mat > 0, 1, mean) * 100,
    Non_Zero_Mean = apply(expr_mat, 1, function(x) mean(x[x > 0]))
  )
  
  # Print summary
  cat("\nMarker Expression Summary:\n")
  print(marker_stats)
  
  # Create distribution plots
  pdf(file.path("output", "stromal_analysis", "plots", "marker_distributions.pdf"))
  par(mfrow=c(3,3))
  for(marker in available_markers) {
    hist(expr_mat[marker,],
         main=paste(marker, "Distribution"),
         xlab="Expression",
         breaks=50)
  }
  dev.off()
  
  return(marker_stats)
}

# Function to identify stromal cells
identify_stromal_cells <- function(se, available_markers) {
  expr_mat <- assay(se, "normcounts")[available_markers,]
  
  # Step 1: VIM+ stromal cells
  vim_positive <- expr_mat["VIM",] > quantile(expr_mat["VIM", expr_mat["VIM",] > 0], 0.25)
  
  # Step 2: PDGFRA/NCAM1 expression (more lenient criteria)
  pdgfra_ncam1_positive <- expr_mat["PDGFRA",] > 0 | expr_mat["NCAM1",] > 0
  
  # Step 3: Low expression of exclusion markers
  negative_markers <- intersect(c("PTPRC", "PECAM1", "EPCAM"), available_markers)
  negative_score <- colMeans(expr_mat[negative_markers,])
  negative_low <- negative_score < quantile(negative_score, 0.75)
  
  # Combine criteria
  stromal_cells <- vim_positive & pdgfra_ncam1_positive & negative_low
  
  # Print summary
  cat("\nCell Population Summary:\n")
  cat("Total cells:", ncol(se), "\n")
  cat("VIM+ cells:", sum(vim_positive), "\n")
  cat("PDGFRA/NCAM1+ cells:", sum(pdgfra_ncam1_positive), "\n")
  cat("Low negative marker cells:", sum(negative_low), "\n")
  cat("Final stromal cells:", sum(stromal_cells), "\n")
  
  return(list(
    is_stromal = stromal_cells,
    vim_positive = vim_positive,
    pdgfra_ncam1_positive = pdgfra_ncam1_positive,
    negative_low = negative_low
  ))
}

# Function to analyze BMP2 expression
analyze_bmp2_expression <- function(se, stromal_results) {
  stromal_cells <- stromal_results$is_stromal
  bmp2_expr <- assay(se, "normcounts")["BMP2", stromal_cells]
  condition <- factor(colData(se)$condition[stromal_cells],
                     levels = c("Control", "ICD"))
  
  # Calculate statistics by condition
  stats <- tapply(bmp2_expr, condition, function(x) {
    c(mean = mean(x),
      median = median(x),
      detection_rate = mean(x > 0) * 100,
      n = length(x))
  })
  
  # Statistical test
  wilcox_test <- wilcox.test(bmp2_expr ~ condition)
  
  # Create visualization
  pdf(file.path("output", "stromal_analysis", "plots", "bmp2_expression.pdf"))
  p1 <- ggplot(data.frame(
    BMP2 = bmp2_expr,
    Condition = condition
  )) +
    geom_violin(aes(x = Condition, y = BMP2, fill = Condition)) +
    geom_boxplot(aes(x = Condition, y = BMP2), width = 0.2) +
    theme_minimal() +
    ggtitle("BMP2 Expression in Stromal Cells")
  print(p1)
  dev.off()
  
  return(list(
    statistics = stats,
    wilcox_test = wilcox_test
  ))
}

# Run Analysis Pipeline ----

# 1. Analyze marker quality
print("Verifying marker availability:")
all_markers <- unlist(markers)
available_markers <- all_markers[all_markers %in% rownames(se)]
missing_markers <- all_markers[!all_markers %in% rownames(se)]
marker_stats <- analyze_marker_quality(se, available_markers)

# 2. Identify stromal cells
stromal_results <- identify_stromal_cells(se, available_markers)

# 3. Analyze BMP2 expression
bmp2_analysis <- analyze_bmp2_expression(se, stromal_results)

# Save results
results <- list(
  marker_stats = marker_stats,
  stromal_results = stromal_results,
  bmp2_analysis = bmp2_analysis
)

saveRDS(results, 
        file.path("output", "stromal_analysis", "processed", 
                  "analysis_results.rds"))

# Additional Visualizations ----

# Create spatial distribution plot
pdf(file.path("output", "stromal_analysis", "plots", "spatial_distribution.pdf"))
ggplot(data.frame(
  x = colData(se)$x_mm,
  y = colData(se)$y_mm,
  Stromal = stromal_results$is_stromal,
  Condition = colData(se)$condition
)) +
  geom_point(aes(x = x, y = y, color = Stromal), size = 0.5) +
  facet_wrap(~Condition) +
  scale_color_manual(values = c("grey80", "red")) +
  theme_minimal() +
  ggtitle("Spatial Distribution of Stromal Cells")
dev.off()

# Create marker correlation plot
pdf(file.path("output", "stromal_analysis", "plots", "marker_correlations.pdf"))
stromal_markers <- c("VIM", "PDGFRA", "NCAM1", "BMP2")
expr_matrix <- t(as.matrix(assay(se, "normcounts")[stromal_markers,]))
cor_matrix <- cor(expr_matrix, method = "spearman")
pheatmap(cor_matrix,
         main = "Marker Correlations",
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

# Results Summary ----
cat("\nFinal Analysis Summary:\n")
cat("=====================\n")

# Print statistical test results
print(bmp2_analysis$statistics)
cat("\nStatistical Test Results:\n")
print(bmp2_analysis$wilcox_test)

# Create summary table
summary_table <- data.frame(
  Metric = c("Total Cells", "Stromal Cells", "BMP2+ Stromal Cells"),
  Control = c(
    sum(colData(se)$condition == "Control"),
    sum(stromal_results$is_stromal & colData(se)$condition == "Control"),
    sum(stromal_results$is_stromal & colData(se)$condition == "Control" & 
        assay(se, "normcounts")["BMP2",] > 0)
  ),
  ICD = c(
    sum(colData(se)$condition == "ICD"),
    sum(stromal_results$is_stromal & colData(se)$condition == "ICD"),
    sum(stromal_results$is_stromal & colData(se)$condition == "ICD" & 
        assay(se, "normcounts")["BMP2",] > 0)
  )
)

print(summary_table)

# Calculate percentages
percentage_table <- data.frame(
  Metric = c("% Stromal of Total", "% BMP2+ of Stromal"),
  Control = c(
    (summary_table$Control[2] / summary_table$Control[1]) * 100,
    (summary_table$Control[3] / summary_table$Control[2]) * 100
  ),
  ICD = c(
    (summary_table$ICD[2] / summary_table$ICD[1]) * 100,
    (summary_table$ICD[3] / summary_table$ICD[2]) * 100
  )
)

cat("\nPercentages:\n")
print(percentage_table)

# Print session information for reproducibility
sessionInfo()