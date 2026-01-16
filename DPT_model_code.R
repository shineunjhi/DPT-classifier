# -------------------------------------------------------
#  Random Forest–based Single-Cell T Cell Classification Model
#
#  This script applies a pre-trained Random Forest classifier
#  to single-cell transcriptomic profiles to estimate class
#  membership probabilities for peripheral CD4⁺, CD8⁺, and
#  CD4⁺CD8⁺ double-positive (DP) T cells.
#
#  The model is trained on transcriptional profiles of
#  sorted memory CD4⁺ and CD8⁺ T cells and enables systematic
#  identification of rare double-positive transcriptional
#  states based on multi-gene expression patterns.
#
#  Model outputs include class-specific probabilities and
#  support downstream evaluation of discrimination and
#  calibration using one-vs-rest strategies.
#
#  Author: Eunji Shin
#  Contact: eunjhi.shin@gmail.com
#  Created: 2025
# -------------------------------------------------------


# -------------------------------------------------------
# Load Data
# -------------------------------------------------------

setwd("/DPT_model")      # Set working directory where data files are stored
load("DPT_model.rda")    # Load pre-trained RF model objects
load("test_matrix.rda")  # Load test dataset


# -------------------------------------------------------
# Automatically Install & Load Required Packages
# -------------------------------------------------------
packages <- c("mlr", "randomForest", "pROC", "caret", "ggplot2", "dplyr")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Load essential libraries
library(mlr)
library(randomForest)
library(pROC)
library(caret)
library(ggplot2)
library(dplyr)



# -------------------------------------------------------
# 1) Prepare Test Matrices
# -------------------------------------------------------
test_matrix$cluster = factor(test_matrix$cluster, levels = c(0,1))
test_classif <- makeClassifTask(data = test_matrix, target = "cluster")


# -------------------------------------------------------
# 2) Compute Random Forest model in Test Set
# -------------------------------------------------------

# Generate predictions for the test dataset using the pre-trained Random Forest model
rfmodel <- predict(rforest, test_classif)

# Extract the underlying Random Forest learner to inspect model details
getLearnerModel(rforest)

# Extract true outcome labels from the prediction object
test_y = as.numeric(as.character(rfmodel$data$truth))

# Extract predicted probabilities for the positive class (class = 1)
RF_test_DP= rfmodel$data$DP

# -------------------------------------------------------
# 3) Performance Evaluation (AUC, ROC)
# -------------------------------------------------------

pred <- rfmodel$data[, c("prob.CD4", "prob.CD8", "prob.DP")]
colnames(pred) <- c("CD4", "CD8", "DP")

multiclass.roc(
  response = rfmodel$data$truth,
  predictor = as.matrix(pred)
)


roc_cd4 <- roc(
  rfmodel$data$truth == "CD4",
  rfmodel$data$prob.CD4
)

roc_cd8 <- roc(
  rfmodel$data$truth == "CD8",
  rfmodel$data$prob.CD8
)

roc_dp <- roc(
  rfmodel$data$truth == "DP",
  rfmodel$data$prob.DP
)

auc(roc_cd4)
auc(roc_cd8)
auc(roc_dp)




# -------------------------------------------------------
# 4) Performance Evaluation (Confusion Matrix)
# -------------------------------------------------------

# Confusion matrix based on RF model
# Sensitivity, specificity, PPV, NPV, F1-score
conf = confusionMatrix(rfmodel$data$response, rfmodel$data$truth, positive = "1")
print(conf$byClass)   # Sensitivity, specificity, etc.
print(conf$overall)   # Accuracy, kappa



# -------------------------------------------------------
# 5) Visualization
# -------------------------------------------------------


RF_df <- data.frame(
  RF = rfmodel$data$prob.DP,
  Outcome = factor(
    ifelse(rfmodel$data$truth == "DP", "DP", "non-DP"),
    levels = c("non-DP", "DP")
  )
)

png("01.DP_RF_prob.png", width = 850, height = 1000)
myplot <-
  ggplot(RF_df, aes(x = Outcome, y = RF, fill = Outcome)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, color = "black") +
  scale_fill_manual(values = c("#4DA3FF", "#FF6B6B")) +  # non-DP vs DP
  labs(
    x = "",
    y = "RF predicted probability (DP)",
    title = "RF predicted probability by DP group"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 30, face = "bold"),
    axis.text.x  = element_text(size = 30, face = "bold"),
    axis.text.y  = element_text(size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, face = "bold")
  )

print(myplot)
dev.off()


# -------------------------------------------------------
# 6) Feature Contribution Plot
# -------------------------------------------------------

# Compute permutation-based feature importance from the Random Forest model
imp_perm = getFeatureImportance(rforest, type = 1)

# Convert feature importance results into a data frame
featureImportance <- data.frame(imp_perm$res)

# Select the top 100 most important features for visualization
featureImportance = featureImportance[order(featureImportance$importance, decreasing = TRUE),]
featureImportance_top100 = featureImportance[1:100,]


coef_plot <- featureImportance_top100 %>%
  arrange(desc(abs(importance))) %>%
  mutate(variable = factor(variable, levels = rev(variable)))  # Preserve order visually

png("02.DPT_features.png", width = 900, height = 1600)
myplot =
  ggplot(coef_plot, aes(x = importance, y = variable, fill = importance > 0)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.85) +
  scale_fill_manual(
    values = c("#FF6B6B", "#4DA3FF"),
    labels = c("Negative effect", "Positive effect")
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    title = "Feature Contribution to RF Model",
    x = "Feature importance score",
    y = "",
    fill = "Effect Direction"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 20, face = "bold"),
    plot.title  = element_text(size = 20, face = "bold"),
    legend.position = "none"
  )
print(myplot)
dev.off()