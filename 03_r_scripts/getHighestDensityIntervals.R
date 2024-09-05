library(tidyverse)
source("../R_notebooks/helperFunctions.R")

results.path <- list.dirs("../01_results/00_panse_fits", recursive = F)
runs <- c("restart_1","final_restart")
files.to.analyze <- crossing(results.path,runs) %>% 
  purrr::pmap(~paste(..1,..2,sep="/"))

files.to.analyze %>% 
  purrr::map(function(x)
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
    
       write_csv(results$Log.NSERate,file=file.path(x,"Parameter_est",nse.files[str_detect(test,"log_scale")]))
       write_csv(results$NSERate,file=file.path(x,"Parameter_est",nse.files[!str_detect(test,"log_scale")]))
  })

