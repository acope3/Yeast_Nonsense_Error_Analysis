Model was run with bash command
```
Rscript --vanilla R_scripts/runPANSEMixModel.R -i 00_data/00_panse_input/2026-03-16_chou_etal_2017_elp1D_1_all_frames_200_ramp.csv -o 01_results/00_panse_fits//2026-03-20_Chou_etal_2017_elp1D_1_all_genes_longer_than_225_codons_200_ramp --dataset Chou_etal_2017_epl1D_1 -d 0 -s 10000 -a 20 -t 5 -n 40 --phi 00_data/00_panse_input/2026-03-16_chou_etal_2017_elp1D_1_all_frames_200_ramp_phi.csv --est_csp --est_hyp --est_phi --normalize_phi --mixture_definition allUnique --max_num_runs 2 --development ~/AnaCoDa_installs/4.4.0/PA_ignore_positions
```
