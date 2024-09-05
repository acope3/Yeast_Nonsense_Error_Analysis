# File descriptions

This directory contains R Notebooks that are used to generate figures and perform post-hoc analyses from model fits of PANSE. These include

- 00_data_quality_and_filtering.Rmd: checking for 5'-ramp and genes that may violate the assumptions of PANSE
- 01_model_comparisons.Rmd: comparing different runs of PANSE to assess how different factors affect parameter estimates e.g., different filtering criteria for which genes are included, inclusion vs. exclusion of first 200 codons on NSE rate estimates
- 02_model_adequacy.Rmd: comparing real and data simulated under posterior parameter estimates, particularly as it relates to the 5'-ramp region
- 03_contextualizing_parameter_estimates.Rmd: comparing PANSE parameter estimates to independent empirical data (e.g., mRNA abundances) and theoretical expectations (e.g., are codons with greater NSE probabilities biased toward the 5'-end)
- 04_testing_for_adaptation_against_nonsense_errors.Rmd: testing for adaptation to reduce the cost of NSEs, including how translation costs (direct, indirect, and NSE costs) vary with gene expression and length