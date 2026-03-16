library(riboWaltz)
library(tidyverse)
library(rtracklayer)

gff_file <- "~/example-datasets/fungi/saccharomyces/annotation/Saccharomyces_cerevisiae_yeast_CDS_w_250utrs.gff3"
bam_file <- "/nobackup/rokaslab/copea1/Public_sequencing/Ribo_seq/Fungi/Scerevisiae/Chou_etal_2017_Mol_Cell/output_trimmed/WT_1/"
length_range <- 28:32

gff <- readGFFAsGRanges(gff_file)

gff_df <- data.frame(transcript=as.character(gff$Name),
                     type = as.character(gff$type),
                     width = width(gff)) %>%
  mutate(type = case_when(
    type == "CDS" ~ "l_cds",
    type == "UTR5" ~ "l_utr5",
    type == "UTR3" ~ "l_utr3",
  )) %>% 
  pivot_wider(id_cols = transcript,names_from=type,values_from=width) %>%
  mutate(l_tr = l_utr5 + l_cds + l_utr3) %>%
  dplyr::select(transcript,l_tr,l_utr5,l_cds,l_utr3) %>%
  as.data.frame()
reads_list <- bamtolist(bamfolder = bam_file, annotation = gff_df)

filtered_list <- length_filter(data = reads_list,
                               length_filter_mode = "custom",
                               length_range = length_range)
psite_offset <- psite(filtered_list, flanking = 6, extremity = "auto")
