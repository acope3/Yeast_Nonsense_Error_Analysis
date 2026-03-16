Model was run with bash command
```
Rscript --vanilla R_scripts/runPANSEMixModel.R -i 00_data/00_panse_input/2025-11-21_ferguson_etal_2023_all_frames_200_ramp.csv -o 01_results/00_panse_fits//2025-11-21_Ferguson_etal_2023_rep1_plus_rep2_all_frames_filter_genes_shorter_than_225_codons_200_ramp_share_nse --dataset Ferguson_etal_2023 -d 0 -s 10000 -a 20 -t 5 -n 48 --phi 00_data/00_panse_input/2025-11-21_ferguson_etal_2023_all_frames_200_ramp_phi.csv --est_csp --est_hyp --est_phi --normalize_phi --share_nse --mixture_definition allUnique --max_num_runs 2 --development ~/AnaCoDa_installs/PA_ignore_positions
```
