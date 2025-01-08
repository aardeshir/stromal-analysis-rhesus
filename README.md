# Rhesus Stromal Cell Analysis

[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](https://choosealicense.com/licenses/mit/)

Analysis of BMP2-expressing stromal cells in primate colon comparing chronic diarrhea vs control cases.

## Project Overview

This project examines BMP2-expressing stromal cells in chronic diarrhea vs control primates, testing parallels with antibiotic-driven pathways. The analysis focuses on comparing samples from:
- ICD case (chronic diarrhea): RB45 (born 2023-04-27)
- Control: RC56.1 (born 2023-05-14)

## Repository Structure
```
├── data/           # Raw and processed data
├── R/              # R scripts for analysis
├── output/         # Analysis outputs
│   ├── plots/      # Generated figures
│   └── processed/  # Processed data files
├── docs/           # Documentation and protocols
└── renv/           # R environment management
```

## Setup

```R
install.packages('renv')
renv::restore()
```

## Usage

Main analysis script: `R/stromal_analysis.R`

## License

[MIT](https://choosealicense.com/licenses/mit/)

## Contact

A. Ardeshir (Ardeshirlab) - aardeshir@tulane.edu