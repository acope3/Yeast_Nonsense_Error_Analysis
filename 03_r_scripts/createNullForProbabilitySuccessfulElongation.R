library(tidyverse)
library(parallel)
library(argparse)
library(AnaCoDa)

parser <- ArgumentParser()
parser$add_argument("-i","--input",help="File storing RFP data",type="character")
parser$add_argument("-o","--output",help="File to output regression parameter estimates",type="character")
parser$add_argument("-n","--num_permutations",help="Number of permutations to perform",type="integer",default=1000)
parser$add_argument("--permutation_type",help="Type of permutation: should be AA, CDS, or Genome to indicate within AA and CDS only, within CDS only, or across CDS",type="character",default="AA")
parser$add_argument("--nserate_file",help="File path to file containing NSE rates from PANSE",type="character")
parser$add_argument("--wait_time_file",help="File path to file containing waiting times from PANSE",type="character")
parser$add_argument("--max_position",help="The maximum position (in terms of codons) to be included in the final regression",type="integer",default=500)
parser$add_argument("--num_cores",help="Number of cores to use to generate mull",type="integer",default=8)
args <- parser$parse_args()
rfp.input <- args$input
output <- args$output
num.perm <- args$num_permutations
permutation.type <- args$permutation_type 
nserate.file <- args$nserate_file
wait.time.file <- args$wait_time_file
max.position <- args$max_position
num.cores <- args$num_cores


rfp <- read_csv(rfp.input)
nserate <- read_csv(nserate.file)
wait.time <- read_csv(wait.time.file)

rfp  <- rfp %>%
  rowwise() %>%
  mutate(AA = codonToAA(Codon)) %>%
  ungroup()


final.df <- mclapply(1:num.perm,function(i){
  if (permutation.type == "AA")
  {
    null <- rfp %>%
      mutate(Codon = with(rfp,ave(Codon,GeneID,AA,FUN=sample)))
  } else if (permutation.type == "CDS"){
    null <- rfp %>%
      mutate(Codon = with(rfp,ave(Codon,GeneID,FUN=sample)))
  } else if (permutation.type == "Genome"){
    null <- rfp %>%
      mutate(Codon = sample(Codon))
  } else {
    stop("Must specify AA, CDS, or Genome for --permutation_type")
  }
  
  null <- null %>%
    left_join(wait.time %>% dplyr::select(Codon,Mean) %>% dplyr::rename(Wait.Time = Mean),by="Codon") %>%
    left_join(nserate %>% dplyr::select(Codon,Mean) %>% dplyr::rename(NSERate = Mean),by="Codon")
  
  total.null <- null %>% 
    mutate(Prob.Successful = (1/NSERate)/(1/NSERate + Wait.Time)) %>%
    filter(Position <= max.position) %>%
    group_by(Position) %>%
    summarize(Total.Prob.Successful = prod(Prob.Successful)^(1/n()))
  
  l <- lm(Total.Prob.Successful ~ log(Position),data=total.null)
  coef.l <- coefficients(l)
  
  df <- data.frame(Permutation=i,Intercept = coef.l[1],Slope=coef.l[2])
  df
}, mc.cores=num.cores
) %>% bind_rows()

write_csv(final.df,output)





