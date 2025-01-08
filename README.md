# Rhesus Stromal Cell Analysis

Analysis of BMP2-expressing stromal cells in primate colon comparing chronic diarrhea vs control cases.

## Project Overview

This project examines BMP2-expressing stromal cells in chronic diarrhea vs control primates, testing parallels with antibiotic-driven pathways. The analysis focuses on comparing samples from:
- ICD case (chronic diarrhea): RB45 (born 2023-04-27)
- Control: RC56.1 (born 2023-05-14)

## Repository Structure

```
├── data/           # Raw and processed data (not tracked in git)
├── R/              # R scripts for analysis
├── output/         # Analysis outputs (not tracked in git)
│   ├── plots/      # Generated figures
│   └── processed/  # Processed data files
├── docs/           # Documentation and protocols
└── renv/           # R environment management
```

## Setup Instructions

1. Clone the repository
2. Install R dependencies:
```R
install.packages('renv')
renv::restore()
```

## Analysis Pipeline

The main analysis script `R/stromal_analysis.R` performs:
1. Data loading and preprocessing
2. Marker quality analysis
3. Stromal cell identification
4. BMP2 expression analysis
5. Statistical testing and visualization

## Contact

A. Ardeshir (Ardeshirlab) - aardeshir@tulane.edu