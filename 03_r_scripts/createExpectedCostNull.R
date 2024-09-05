library(tidyverse)
library(ggpubr)
library(feather)
library(AnaCoDa)


suppressMessages(source("../R_notebooks/helperFunctions.R"))


deta <- read_csv("../Data/ROC_Estimates/saccharomyces_cerevisiae.max.cds_Selection.csv") %>% 
  dplyr::rename(DEta = Mean)
rfp <- read_csv("../Data/PANSE_Formatted/Scerevisiae/Weinberg_etal_2016/weinberg_etal_2016_all_frames.csv")
if ("gene" %in% colnames(rfp))
{
  rfp <- rfp %>% dplyr::rename(GeneID=gene)
}
rfp<- rfp %>%
  rowwise() %>%
  mutate(AA=codonToAA(Codon)) %>%
  ungroup()
permutations <- 1000
results.files <- c("../Results/PANSE/Scerevisiae/Weinberg/All_frames/Full_genome/Ramp/2023-07-20_Weinberg_etal_all_genes_filter_genes_shorter_than_225_codons_200_ramp_run_2/final_restart/")


init.cost <- c(4)
a.2 <- 4
w.0 <- 10.8
Ne <- 13600000 ## S. paradoxus effective pop. size
q <- 4.19*10^-7 ## from Gilchrist 2007 MBE, but this does include NSE cost
target.mean.elong.rate <- 9.3 #codon/s based on Shah et al. 2013 and Arava et al. 2013, Dao-Duc and Song 2018
deta.file <- "../Data/ROC_Estimates/saccharomyces_cerevisiae.max.cds_Selection.csv"
index <- 1

panse.result <- getParameterDataFrames(results.files[1])
wait.time.panse <- panse.result$Wait.Time %>%
  dplyr::rename(Wait.Time = Mean)

wait.time.unit.codon.per.sec <- rescalePANSEWaitingTimes(wait.time.panse,
                                                         target.mean.elong.rate = target.mean.elong.rate)

dist.C <- calculateCDistribution(wait.time.unit.codon.per.sec,
                                 Ne = Ne,
                                 q = q,
                                 deta.file = deta.file)


C.roc <- calculateMedianC(dist.C,include.discrepancy.codons = T)
C.assembly <- 5.5
C.vec <- c(C.roc,C.assembly)

df.list <- vector(mode="list",length = length(results.files) * length(init.cost) * length(C.vec))

for (result in results.files)
{
  for (a.1 in init.cost)
  {
    for (C in C.vec)
    {
      rfp.perm.all <- purrr::map(1:permutations,function(x){
        rfp.perm <- rfp %>% 
          mutate(Codon = with(rfp,ave(Codon,GeneID,AA,FUN=sample)))
        rfp.perm.exp.cost <- expectedCost(panse.result.file.path = result,
                                          rfp = rfp.perm,
                                          a.1 = a.1,
                                          a.2 = a.2,
                                          w.0 = w.0,
                                          C = C,
                                          target.mean.elong.rate = target.mean.elong.rate)
        rfp.perm.exp.cost <- rfp.perm.exp.cost %>%
        dplyr::select(c(GeneID,
                         Expected.cost.overall,
                         Result.file.path,
                         Initiation.cost.a1,
                         Elongation.cost.a2,
                         w.0,
                         C,
                         Target.elongation.rate)
         )
        return(rfp.perm.exp.cost)
      }) %>% 
        purrr::reduce(left_join,by=c("GeneID","Result.file.path","Initiation.cost.a1","Elongation.cost.a2","w.0","C","Target.elongation.rate"))
      df.list[[index]] <- rfp.perm.all
      index <- index + 1
    }
  }
}

df <- df.list %>% bind_rows()
write_feather(df,"../Paper/Data/2024-07-17_expected_nse_cost_null_distributions_cost_per_nse_and_overall_2023_07_20_post_200_fit_full_gene_with_proper_waiting_time_units.feather")


