Model was run with bash command
```
Rscript --vanilla R_scripts/runPANSEMixModel.R -i 00_data/00_panse_input/2025-11-05_wu_etal_2019_all_frames_21nt_28nt_200_ramp.csv -o 01_results/00_panse_fits//2025-11-05_Wu_etal_2019_28nt_21nt_rfp_all_frames_filter_genes_shorter_than_225_codons_200_ramp --dataset Wu_etal_2019 -d 0 -s 10000 -a 20 -t 5 -n 48 --phi 00_data/00_panse_input/2025-11-05_wu_etal_2019_all_frames_21nt_28nt_phi.csv --est_csp --est_hyp --est_phi --normalize_phi --mixture_definition allUnique --max_num_runs 2 --development ~/AnaCoDa_installs/PA_ignore_positions
```
