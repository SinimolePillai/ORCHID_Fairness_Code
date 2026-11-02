# Pipeline-Aware Fairness Auditing for Sequential Healthcare Decision Systems: Evidence from the Organ Donation Process

## Overview

This repository contains the R analysis code supporting the manuscript *“Pipeline-Aware Fairness Auditing for Sequential Healthcare Decision Systems: Evidence from the Organ Donation Process”*, submitted to the International Journal of Medical Informatics.

The study examines predictive performance, calibration and demographic disparities across three sequential stages of the organ donation process:

1. Approach
2. Authorisation
3. Procurement

The analysis uses logistic regression (LR) and XGBoost to evaluate stage-specific predictions and estimated cumulative pipeline reachability. It also examines how differences across stages may accumulate over the decision pathway.

The study is observational. Estimated demographic disparities should not automatically be interpreted as evidence of discrimination.

## Data Source

The analysis uses the Organ Retrieval and Collection of Health Information for Donation (ORCHID) dataset, available through PhysioNet.

The dataset is not included in this repository. Researchers must obtain the data independently and comply with the applicable PhysioNet access and data-use requirements.

## Software Requirements

The analysis was developed using R version 4.5.1.

Required R packages:

- `data.table`
- `dplyr`
- `ggplot2`
- `pROC`
- `xgboost`
- `Matrix`

Package versions should be checked against the accompanying session information where available.

## Repository Structure

**`ORCHID_code.R`**

Consolidated R analysis script containing data preprocessing, cohort construction, stage-specific LR and XGBoost modelling, calibration, predictive performance evaluation, subgroup fairness assessment, cumulative pipeline reachability estimation, bootstrap uncertainty analysis and sensitivity analyses.

**`Sessioninfo.txt`**

R environment and package information. 

**`LICENSE`**

Repository licence information.

## Analytical Specifications

The primary analysis uses a temporal development, calibration and test design:

- Development: 2015–2019
- Calibration: 2020
- Test: 2021

The primary reduced-predictor specification includes age, gender, race and organ procurement organisation (OPO).

Additional analyses examine expanded predictor specifications, complete-case restrictions and random data partitions.

Logistic regression is fitted using maximum likelihood. XGBoost predictions are calibrated using Platt scaling.

The analysis includes stage-specific predictive performance, calibration, subgroup disparities and estimated cumulative reachability. Bootstrap confidence intervals are conditional on the fitted models and do not incorporate model-refitting uncertainty.

## Reproducibility

To reproduce the analysis:

1. Obtain authorised access to the ORCHID dataset through PhysioNet.
2. Download the relevant ORCHID referral data.
3. Install the required R packages.
4. Download `ORCHID_code.R`.
5. Configure the input data location and output directory as prompted or specified in the script.
6. Execute the script in R.

The script generates analytical outputs for the primary and sensitivity specifications.

## Interpretation and Limitations

The analysis is based on observational referral data. Differences in predicted outcomes or estimated cumulative reachability do not, by themselves, establish unjustified discrimination.

Clinical eligibility and other potentially relevant factors may not be fully represented. The estimated cumulative reachability measure is a modelling approximation and should not be interpreted as an observed individual-level probability of procurement.

The repository does not provide validated causal estimates or an experimentally demonstrated pipeline-aware fairness intervention.

## Data Availability

The ORCHID dataset is available through PhysioNet, subject to its applicable access conditions. Individual-level data are not redistributed in this repository.

## Contact

Dr. Sinimole K. R.  
Rajagiri Business School, Kochi, India
