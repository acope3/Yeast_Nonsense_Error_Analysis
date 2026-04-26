Model was run with bash command
```
Rscript --vanilla R_scripts/runPANSEMixModel.R -i 00_data/00_panse_input/2026-03-24_weinberg_etal_2016_replicate_all_frames_w_3_utrs_200_ramp.csv -o 01_results/00_panse_fits//2026-03-25_Weinberg_etal_2016_all_genes_stop_codons_3_utrs_filter_genes_shorter_than_225_codons_200_ramp_0.2_upper_limit --dataset Weinberg_etal_2016_3_utrs -d 0 -s 10000 -a 20 -t 5 -n 37 --phi 00_data/00_panse_input/2025-10-07_weinberg_etal_2016_all_frames_phi.csv --est_csp --est_hyp --est_phi --normalize_phi --include_stop --nserate_uniform_upper_limit 0.2 --mixture_definition allUnique --max_num_runs 2 --development ~/AnaCoDa_installs/Allow_stop_codons
```
