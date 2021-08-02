# An entropy-based approach to measuring morphological productivity

Data and code for Elizabeth Pankratz's 2021 Universität Potsdam MSc thesis "Morphological productivity: A Bayesian modelling approach" under the supervision of Shravan Vasishth and Titus von der Malsburg.

## Contents

The structure of this repository largely mirrors the structure of the thesis.
Each directory (except `1_data/`) contains the code used for the analyses in the section with the corresponding number.

- `1_data/`: Create the samples used throughout the thesis.
- `2_interpretability/`: Bootstrap the three large samples that were used to illustrate entropy and alternate productivity measures.
- `3_validity/`: Identify German simplexes, compute the three productivity factors used in the Bayesian linear regression analysis, conduct analysis.
- `4_applicability/`: Bootstrap the 35 smaller samples and the samples from RIDGES and evaluate stability of entropy scores.
- `5_outlook/`: Find a mathematical law that approximates entropy as a function of bootstrapped sample size for the three large samples.

Each directory contains further READMEs that detail its contents.
