library(tidyverse)
source("../02_r_notebooks//helperFunctions.R")

results.path <- list.dirs("../01_results/00_panse_fits/2026-01-07_Weinberg_etal_2016_all_frames_filter_genes_shorter_than_225_codons_200_ramp_updated_posterior_variable_nse/", recursive = F)
#runs <- c("restart_1","final_restart")
#files.to.analyze <- crossing(results.path,runs) %>% 
#  purrr::pmap(~paste(..1,..2,sep="/"))
files.to.analyze <- results.path

files.to.analyze %>% 
  purrr::map(function(x)
    {
       if (file.exists(file.path(x,"R_objects","parameter.Rda")))
       {
         results <- getParameterDataFrames(x)
         nse.hdi <- getHDI(file.path(x,"R_objects","parameter.Rda"),log = F)
         log.nse.hdi <- getHDI(file.path(x,"R_objects","parameter.Rda"),log = T)
         results$Log.NSERate <- results$Log.NSERate %>%
           left_join(log.nse.hdi,by="Codon") %>%
           mutate(`2.5%` = HDI_low,
                  `97.5%` = HDI_high) %>%
           dplyr::select(-HDI_low,-HDI_high, -contains("Median"))
         results$NSERate <- results$NSERate %>%
           left_join(nse.hdi,by="Codon") %>%
           mutate(`2.5%` = HDI_low,
                  `97.5%` = HDI_high) %>%
           dplyr::select(-HDI_low,-HDI_high,-contains("Median"))
         nse.files <- list.files(file.path(x,"Parameter_est"),pattern="NSERate")
      
         write_csv(results$Log.NSERate,file=file.path(x,"Parameter_est",nse.files[str_detect(nse.files,"log_scale")]))
         write_csv(results$NSERate,file=file.path(x,"Parameter_est",nse.files[!str_detect(nse.files,"log_scale")]))
       }
  })

