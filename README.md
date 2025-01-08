# Stromal Cell Analysis

Analysis of BMP2-expressing stromal cells in primate colon comparing chronic diarrhea vs control cases.

## Overview

This project examines BMP2-expressing stromal cells in chronic diarrhea vs control primates:
- ICD case (chronic diarrhea): RB45 (born 2023-04-27)
- Control: RC56.1 (born 2023-05-14)

## Requirements

- R 4.4.0
- SummarizedExperiment 1.34.0
- SpatialExperiment 1.14.0
- scater 1.32.1
- ggplot2 3.5.1
- tidyverse packages

## Usage

Main analysis script: `R/stromal_analysis.R`

## Data

Place required files in `data/` directory:
- Additional_metadata.csv
- PR03_Chronic/data/processed/:
  - counts.RDS
  - xy.RDS
  - negcounts.RDS 
  - metadata.RDS

## License

[MIT](LICENSE)

## Contact

A. Ardeshir (ArdeshirLab) - aardeshir@tulane.edu