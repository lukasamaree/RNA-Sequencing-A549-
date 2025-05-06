
# Lung Carcinoma Single-Cell RNA-Seq Analysis

This project applies unsupervised machine learning and statistical visualization techniques to uncover cellular heterogeneity and intercellular communication in lung carcinoma using single-cell RNA sequencing data.

## Overview

Single-cell RNA-seq (scRNA-seq) data allows for high-resolution analysis of gene expression at the cellular level. In this project, we:

- **Cleaned and visualized** over **1,000 RNA sequences** using violin plots, Principal Component Analysis (PCA), and K-Means clustering.
- **Identified 7 distinct cell types** involved in lung carcinoma through unsupervised clustering techniques.
- **Mapped intercellular communication pathways** to better understand the tumor microenvironment and signaling interactions among cancer cell subtypes.

## Key Methods and Tools

- **Language**: R  
- **Data Visualization**: Violin plots, PCA plots, cluster heatmaps  
- **Dimensionality Reduction**: Principal Component Analysis (PCA)  
- **Clustering**: K-Means Clustering  
- **Biological Insights**: Tumor heterogeneity, cell differentiation, intercellular signaling

## Findings

- Revealed **7 heterogeneous cell clusters** that likely represent unique biological cell types within the tumor.
- Showed **variation in gene expression patterns** across cell types, aiding in classification and annotation.
- Identified **intercellular signaling networks**, shedding light on possible therapeutic targets in the tumor microenvironment.

## Impact

This analysis supports improved interpretation of lung cancer at the single-cell level and contributes to potential advancements in **precision oncology**, particularly in understanding **cellular diversity and communication within tumors**.

## Future Work

- Integrate ligand-receptor interaction databases for deeper signaling analysis.
- Use graph-based clustering (e.g., Louvain or Leiden) for improved cluster resolution.
- Apply pseudotime trajectory analysis to track cell lineage differentiation.


# Images

<img width="942" alt="image" src="https://github.com/user-attachments/assets/f7a9d3eb-9ee8-46f7-bfc6-07fb92a6bb59" />

I used PCA, SNN, and UMAP to reduce the dimension of the RNA sequences and plot group and plot them. I then used most the active genes in each cluster to label them by the cells in which those genes are most prominent.

<img width="463" alt="image" src="https://github.com/user-attachments/assets/26b3095c-1fc2-4405-b8b4-150db2d05d0e" />

Using the clustered data I used average gene expression in each cluster to identify which clusters interacted with each other. This is due to the fact if expression levels are high their is a higher probability that it interacted with other cells in comparison to cells that do not have high gene expression.


