# 3D-Gene-Clustering
### Triclustering model for three-dimensional time-series gene expression data

## Introduction
This repository provides the R implementation of a three-dimensional gene expression data clustering technique. The model is based on a multivariate Gaussian mixture model (GMM) under a maximum likelihood framework, integrated with Legendre polynomials to simulate three-dimensional patterns and the SAD model to simulate the covariance structure. 

This tool is designed to identify gene clusters that share similar expression patterns across three dimensions: Time, Tissue (Space), and Environment.

## Repository Structure
- funs.R: Core function library containing the mathematical implementation of the GMM and triclustering algorithms.
- allfw.R: Main execution script to run the clustering process.
- quan.R: Data analysis and hypothesis testing scripts used in the study.
- data/: Processed datasets, including simulated data (cluster_data.txt) and empirical Arabidopsis thaliana data.

## Requirements
To run the scripts, you need R (version >= 4.0.0) and the following packages:
- MASS
- mvtnorm

## Usage
1. Clone the repository:
   git clone https://github.com/DBYLqk/3D-Gene-Clustering.git

2. Run the analysis:
   Open R or RStudio, set the working directory to the repository folder, and source the main script:
   source("funs.R")
   source("allfw.R")
   
   Note: Ensure the data file paths in the scripts match the location of files in the "data/" folder.

## Citation
If you use this code or model in your research, please cite our paper:
QianKun Liu, MengYuan Zhu, Dongchao Ji, Libo Jiang. Triclustering model for three-dimensional time-series gene expression data. (2026), International Journal of Molecular Sciences.

## License
This project is licensed under the MIT License - see the LICENSE file for details.
