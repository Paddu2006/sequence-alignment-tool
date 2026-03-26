# Sequence Alignment Tool

A Python implementation of two fundamental bioinformatics alignment
algorithms with BLAST-like output and visualization.
Built as Phase 2 - Project 1 of my bioinformatics-to-tech portfolio.

## Algorithms Implemented
- Needleman-Wunsch — Global sequence alignment
- Smith-Waterman — Local sequence alignment

## Features
- Global alignment using Needleman-Wunsch algorithm
- Local alignment using Smith-Waterman algorithm
- BLAST-like output format
- Scoring matrix heatmap visualization
- Alignment statistics bar chart
- Supports multiple sequence pairs
- Auto-saves results to timestamped file

## Scoring Scheme
- Match    : +2
- Mismatch : -1
- Gap      : -2

## Sample Output
Alignment Type : Global (Needleman-Wunsch)
Score          : 26
Identities     : 14/16 (87.5%)
Mismatches     : 2/16
Gaps           : 0/16

Query  1  ATGCTAGCTAGCATGC  16
          ||||.|||||||||.|
Sbjct  1  ATGCAAGCTAGCATCC  16

## Tech Stack
- Python 3.13
- NumPy
- Matplotlib
- VS Code
- Git + GitHub

## Usage
python alignment.py

Enter two DNA sequences and choose alignment type:
1. Global (Needleman-Wunsch)
2. Local (Smith-Waterman)
3. Both

## Project Roadmap
- Phase 2 Project 1 - Sequence Alignment Tool - Done
- Phase 2 Project 2 - Gene Expression Analysis - Coming soon
- Phase 2 Project 3 - Protein Structure Analyzer - Coming soon
- Phase 2 Project 4 - Genomic Variant Database - Coming soon
- Phase 2 Project 5 - Disease Prediction ML - Coming soon

## Author
Padma Shree Jena
Bioinformatics + Tech Enthusiast | Python | R | Bash
GitHub: https://github.com/Paddu2006

## License
MIT License
