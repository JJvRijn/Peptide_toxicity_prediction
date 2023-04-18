# Pipeline feature_calculation

> *Note:* This is a `README.md` boilerplate generated using `Kedro 0.18.7`.

## Overview

This pipeline first grabs the unique sequences from the dataset, so that no feautures get calculated twice for the same peptide.
Then the pipeline calculates different kind of features for the sequences and stores them with the sequence attached for easy merger.


## Pipeline inputs

The complete raw data
But actually any dataframe with a row called seq which is filled with strings containing aminoacid-sequences

## Pipeline outputs

5 dataframes containing different kind of features
