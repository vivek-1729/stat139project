# Stat 139 Project

Al Liu, Eric Jiang, Leonard Mauvernay, Vivek Shah

We study whether the slope linking treatment effects on overall survival to those on progression-free survival differs for immune checkpoint inhibitors. The data come from randomized phase III trials in solid tumors.

## Data

- data/raw/ctgov holds the trial results from ClinicalTrials.gov.
- data/raw/hemonc_2026-08-28 holds HemOncKB version 54 (Harvard Dataverse, doi:10.7910/DVN/FPO4HB, license CC BY-NC-SA 4.0). We use it for the drug classes.
- data/processed holds the files our code builds.

## Code

- code/build_data.R builds the analysis dataset from the raw files.
- eda/eda.Rmd is the EDA milestone. It sources the script.

Knit eda.Rmd from the eda folder or run Rscript build_data.R from the code folder

eda.Rmd also needs the R packages sandwich and lmtest.
