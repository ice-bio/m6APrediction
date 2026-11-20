# m6APrediction

Random Forest based R package for m6A RNA modification site prediction

## Overview

`m6APrediction` is an R package built for predicting m6A RNA modification sites using a pre-trained random forest classifier. Two core functions are provided: 1. `prediction_multiple()`: Batch prediction from a feature table containing multiple candidate sites 2. `prediction_single()`: Prediction for one individual RNA candidate site

The model was trained on RNA sequence and genomic features.

## Installation

Install directly from GitHub with `remotes`:

```{r}
# install remotes if needed
install.packages("remotes")
remotes::install_github("ice-bio/m6APrediction")

library(m6APrediction)

# Load the pre-trained model
rf_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))

### Batch prediction
input_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
batch_result <- prediction_multiple(rf_fit, input_df, positive_threshold = 0.5)
head(batch_result)

### Single site prediction
single_result <- prediction_single(
  ml_fit = rf_fit,
  gc_content = 0.6,
  RNA_type = "mRNA",
  RNA_region = "CDS",
  exon_length = 12,
  distance_to_splice = 5,
  evolutionary_conservation = 0.8,
  DNA_5mer = "GGACA",
  positive_threshold = 0.5
)
single_result
```

## Output explanation

pred_prob: predicted probability of being an m6A modified site pred_label: classification result, "Positive" = m6A site, "Negative" = non-m6A site
