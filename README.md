# Cancer Prevention Knowledge and Adherence to Lifestyle Recommendations: A Longitudinal Study

## Overview

This repository contains the R analysis code for a longitudinal study examining associations between cancer prevention knowledge and adherence to lifestyle recommendations among adults in Chongqing, China.

The analyses include descriptive statistics, cross-lagged panel models (CLPM), sensitivity analyses, and Gaussian graphical network analyses.

## Code Structure

The analysis consists of five scripts:

- `00_data_cleaning.R` — Data preparation, scoring, and analytic sample selection.
- `01_descriptive_analysis.R` — Descriptive statistics, longitudinal changes, and inverse probability weighting.
- `02_clpm_analysis.R` — Primary and domain-specific cross-lagged panel models.
- `03_sensitivity_analysis.R` — Alternative adherence scores, follow-up restrictions, IPW models, and subgroup analyses.
- `04_network_analysis.R` — Network estimation, community detection, bridge centrality, bootstrap stability, and network comparison tests.

The `checks/` directory contains supplementary scripts for verifying selected statistical outputs. The `environment/` directory contains information about the R computing environment.

## Software Requirements

Analyses were conducted using R 4.6.1.

Required R packages include:

`dplyr`, `haven`, `lubridate`, `gtsummary`, `flextable`, `officer`, `ggplot2`, `lavaan`, `ragg`, `moments`, `survey`, `cobalt`, `boot`, `bootnet`, `qgraph`, `networktools`, `NetworkComparisonTest`, `igraph`, `tidyr`, `patchwork`, `gridExtra`, and `scales`.

## Running the Analyses

From the project directory, run:

`Rscript run_all.R`

Alternatively, execute scripts `00` through `04` sequentially in RStudio.

Input files should be placed in the `data/` directory. The input directory can also be specified using the `CEBP_DATA_DIR` environment variable.

Results are saved in the `output/` directory by default. An alternative output directory can be specified using `CEBP_OUTPUT_DIR`.

## Data Availability and Reproducibility

The participant-level SAS datasets (`base.sas7bdat` and `follow.sas7bdat`) are not publicly available because of privacy and institutional data-sharing restrictions.

The R scripts require datasets prepared through an upstream SAS preprocessing workflow. The SAS preprocessing scripts are not included in this repository and may be requested from the authors, subject to institutional approval.

The analysis code is provided to document the statistical procedures and facilitate methodological assessment and reuse. Reproducing the published numerical results requires access to the restricted study datasets.