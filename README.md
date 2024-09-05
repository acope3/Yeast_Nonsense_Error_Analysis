# Yeast_Nonsense_Error_Analysis
Analysis of nonsense errors using yeast ribosome profiling data from Weinberg et al. Cell Reports 2016. 

## Ribosome profiling
 
Ribosome profiling data from Weinberg et al. 2016 Cell Reports (SRA SRR1049521) was downloaded and analyzed using the [riboviz 2](https://github.com/riboviz/riboviz/) pipeline. 
The relevant files for recreating this analysis can be found in [example-datasets](https://github.com/riboviz/example-datasets/tree/main/fungi/saccharomyces).
Specifically, riboviz 2 was run using the YAML configuration file `Weinberg_2016_RPF_1_sample_cerevisiae_CDS_w_250utrs_config_2-1.yaml`. 
See the riboviz 2 documentation for instructions on installing and running riboviz 2.

## Fitting the model

PANSE was fit using a development version of [AnaCoDa](https://github.com/acope3/RibModelFramework/tree/develop). 
For commands used to perform each analysis in `01_results/00_panse_fits/`, a `README.md` is included with the bash command used to run AnaCoDa. 
File paths will need to be altered.
Do not use the CRAN version, which has yet to be updated to reflect the changes made to PANSE.
Please note that PANSE is slow, often taking days to complete (see mcmc_runtime.csv files for each run), even with a relatively large number of threads (48-72 typically used).
Do not run on your desktop. 

## Post model fit analyses

To recreate the analyses and figures in the paper, see `02_r_notebooks` which contains multiple R Markdown files.


