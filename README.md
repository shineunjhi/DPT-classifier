# DPT-classifier
The DPT classifier is a supervised machine learning model trained on scRNA-seq data to accurately classify CD4⁺, CD8⁺, and rare DPT cells. Using a random forest algorithm with optimized parameters, the model achieved 99% accuracy, enabling robust identification of T cell subsets at single-cell resolution.


# Model Design and Training
The model is based on a random forest algorithm and was trained using scRNA-seq data obtained from MACS-sorted memory T cells. To overcome the class imbalance caused by the low abundance of DPT cells, the training process incorporated undersampling of the more frequent CD4⁺ and CD8⁺ populations, ensuring balanced representation across all three T cell subsets.

A total of 5,000 hyperparameter combinations were systematically explored using grid search and 3-fold cross-validation, optimizing the model’s accuracy, precision, recall, and F1-score. Parameters such as the number of decision trees, maximum depth, and feature selection criteria were fine-tuned to achieve optimal performance.

# Performance and Accuracy
The final model demonstrated exceptional classification performance:
1) 99% overall accuracy across CD4⁺, CD8⁺, and DPT T cell subsets
2) 96.5% accuracy specifically for identifying DPT cells

This level of accuracy confirms the model’s ability to detect rare DPT cells with high fidelity, while maintaining strong performance for more abundant T cell types.

# Biological Validation and Application
Predicted DPT cells were validated by their unique transcriptomic and TCR profiles, including:
1) Concurrent expression of CD4 and CD8 genes
2) Enrichment for cytotoxic molecules (e.g., GZMB, GNLY, NKG7)
3) Elevated levels of chemokine genes (e.g., CCL5, CX3CR1)
4) Increased incidence of identical TCR clonotypes and clonal expansion, suggestive of antigen specificity and memory-like features

The model was also successfully applied to external single-cell datasets, including COVID-19 patient samples, achieving 95% classification accuracy and revealing dynamic roles for DPT cells during infection and recovery phases. These findings highlight the model's robustness and applicability across biological conditions and independent datasets.
