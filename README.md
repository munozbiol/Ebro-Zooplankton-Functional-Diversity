# Zooplankton Functional Diversity in Ebro Reservoirs

This repository contains the R statistical workflow to analyze **zooplankton functional diversity** in the Ebro River basin reservoirs. The primary objective is to evaluate how community functional traits respond to environmental, spatial, and temporal gradients.

## 📋 Project Summary

The study utilizes a trait-based ecology approach to go beyond classical taxonomy. It analyzes abundance data and functional traits (feeding type, body weight, reproduction form, etc.) alongside complex physicochemical variables.

## 🛠️ Workflow

The R script is organized into the following critical stages:

### 1. Data Preprocessing and Curation
*   **Taxonomy Cleaning:** Standardizing names between abundance and trait matrices to ensure compatibility with the `FD` package.
*   **Data Imputation:** Using Random Forest (`missForest`) to handle missing values in environmental variables robustly.
*   **Normalization:** Logarithmic transformation of environmental variables and collinearity correction using **VIF (Variance Inflation Factor)**.

### 2. Diversity Indices Calculation
*   **Functional Diversity (Alpha):** Calculation of Functional Richness (FRic), Evenness (FEve), Divergence (FDiv), Dispersion (FDis), and Rao's Entropy (RaoQ).
*   **Taxonomic Diversity:** Shannon, Simpson, and Species Richness indices via `vegan`.

### 3. Multivariate Analysis
*   **PCA (Principal Component Analysis):** Visualization of main environmental gradients (eutrophication and mineralization) grouped by trophic state and reservoir type.
*   **RDA (Redundancy Analysis):** Direct relationship between functional groups and environmental variables using **Hellinger** transformation.

### 4. Variable Selection (Machine Learning)
*   Use of **Conditional Random Forest (`cforest`)** and conditional permutation importance (`permimp`) to identify the environmental predictors with the greatest influence on each diversity index, avoiding overfitting.

### 5. Advanced Ecological Modeling (GAMMs)
*   Implementation of **Generalized Additive Mixed Models (GAMMs)** to capture non-linear relationships.
*   **Fixed Effects:** Selected environmental variables (Conductivity, Chlorophyll-a, Suspended Solids, etc.).
*   **Random Effects:** Modeling the hierarchical structure of the data by including *Reservoir*, *Year*, *Location*, and *Type* as random factors to control for spatial and temporal dependence.

## 📦 Key Libraries

The analysis relies on the following R libraries:
*   `FD`: Calculation of functional indices.
*   `vegan`: Community ecology and ordinations.
*   `mgcv` & `itsadug`: Non-linear modeling (GAMMs).
*   `party` & `permimp`: Random Forest for variable selection.
*   `ggplot2`, `factoextra` & `ggrepel`: Advanced data visualization.

## 📊 Key Findings

The analysis reveals that:
1.  **Conductivity** and **trophic state** are the primary drivers of functional structure.
2.  There is significant **interannual** variability, suggesting the influence of stochastic events (droughts/floods).
3.  Reservoir morphometry (depth and volume) plays a secondary but relevant role in functional richness.

---
*This project is part of the research for the AIL congress and associated publications at the Universitat de València.*

**Contact:**
If you have questions regarding the implementation or the data, please open an *issue* in this repository.
