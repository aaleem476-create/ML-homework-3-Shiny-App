# TB NanoString Shiny App

## Description
This Shiny application explores clustering patterns in the TB NanoString dataset using three dimensionality reduction methods:

- Principal Component Analysis (PCA)
- Uniform Manifold Approximation and Projection (UMAP)
- Non-negative Matrix Factorization (NMF)

Users can select a dimensionality reduction method and choose the number of clusters (k). The app displays the resulting clusters, cluster centers, a cross-tabulation with TB status, and the Adjusted Rand Index (ARI).

## Files
- `app.R` — Shiny application code.
- `TBnanostring.rds` — Dataset used by the application.

## Requirements
Install the required R packages:

```r
install.packages(c("shiny", "ggplot2", "dplyr", "uwot", "NMF", "mclust"))
```

## How to Run

1. Download or clone this repository.
2. Open the project folder in RStudio.
3. Make sure `app.R` and `TBnanostring.rds` are in the same folder.
4. Open `app.R`.
5. Click **Run App** in RStudio.

Alternatively, run:

```r
shiny::runApp()
```