# gdT Classification from scRNA-seq

Transcriptome-based prediction of γδ T cells (gdT) from single-cell RNA sequencing data without relying on TCR sequencing.

This project presents a simple and interpretable framework using Principal Component Analysis (PCA) and Logistic Regression to identify rare T cell populations.


---

## Overview

γδ T cells (gdT) are a rare and heterogeneous T cell population that are difficult to identify using transcriptomic data alone.

Most existing approaches rely on TCR sequencing, which is not always available in large-scale single-cell datasets.

This project proposes a simple and interpretable alternative using:
- Principal Component Analysis (PCA)
- Logistic Regression

to classify gdT versus αβ T cells (abT) using only gene expression data.


---

## Datasets

- **Training dataset:** Song et al. (2023), PBMC scRNA-seq dataset with paired TCR annotations  
  https://zenodo.org/records/7989561

- **External validation dataset:** Domínguez Conde et al. (2022) tissue-atlas dataset  
  https://www.tissueimmunecellatlas.org/

Please download the datasets manually and place them in the `data/` directory before running the analysis.


---

## Methodology

1. Data preprocessing (normalization and log transformation)  
2. Feature selection using Highly Variable Genes (HVGs)  
3. Dimensionality reduction via PCA (computed on HVGs)  
4. Classification using Logistic Regression  
5. Evaluation using ROC AUC and Precision-Recall metrics  


---

## Results

- ROC AUC: ~0.97  
- PR AUC: ~0.58 (baseline ~0.009)

The model demonstrates strong ranking performance (ROC AUC), while maintaining meaningful precision-recall performance under severe class imbalance.

These results indicate that transcriptomic signals alone can partially recover gdT identity, even without TCR information.


---

## Key Insights

- Simple, interpretable models can achieve competitive performance compared to more complex approaches  
- Precision-Recall metrics are essential for evaluating rare cell populations  
- Gene expression alone contains sufficient signal to identify gdT cells to a meaningful extent  


---

## Author

Jaeyoung Shim
BbiomedSc, Biomedical Sciences  
The University of Hong Kong (HKU)
