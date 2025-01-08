# Analysis Protocol

## Data Processing Pipeline

### 1. Quality Control
- Cell filtering based on negative control probes
- Marker validation using expression thresholds

### 2. Stromal Cell Identification
- VIM+ cell selection
- PDGFRA/NCAM1 expression validation
- Exclusion of contaminating lineages

### 3. BMP2 Expression Analysis
- Normalization method: Library size normalization
- Statistical testing: Wilcoxon rank-sum test
- Visualization: Violin plots with overlaid boxplots

## Reproducibility
- R version 4.4.0
- All package versions specified in renv.lock
- Random seed set for reproducible results