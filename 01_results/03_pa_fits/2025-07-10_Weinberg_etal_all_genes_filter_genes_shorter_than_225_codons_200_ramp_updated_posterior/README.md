Model was run with bash command
```
Rscript --vanilla R_scripts/runPANSEMixModel.R -i 00_data/00_panse_input/2023-06-12_weinberg_etal_2016_all_frames_200_ramp.csv -o 01_results/03_pa_fits//2025-07-10_Weinberg_etal_all_genes_filter_genes_shorter_than_225_codons_200_ramp_updated_posterior --dataset Weinberg_etal_2016 -d 0 -s 10000 -a 20 -t 5 -n 48 --restart_file 01_results/03_pa_fits/2025-07-10_Weinberg_etal_all_genes_filter_genes_shorter_than_225_codons_200_ramp_updated_posterior/restart_1/Restart_files/rstartFile.rst_final --est_csp --est_hyp --est_phi --ignore_nse --mixture_definition allUnique --max_num_runs 2 --development ~/AnaCoDa_installs/4.4.0/PA_ignore_positions
```
