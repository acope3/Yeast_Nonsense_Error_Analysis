library(tidyverse)
library(rhdf5)
library(Biostrings)
library(rtracklayer)

source(file.path("/home","copea1","riboviz","rscripts", "read_count_functions.R"))
source(file.path("/home","copea1","riboviz","rscripts", "stats_figs_block_functions.R"))

#' CalcAsiteFixed(): Calculate read A-site using a fixed displacement for fixed read lengths
#'
#' The assignment rules are specified in a user-supplied data frame, `asite_displacement_length`.
#'
#' @param reads_pos_length matrix of read lengths and positions (e.g. as given by GetGeneDatamatrix(gene, dataset, hd_file) )
#' @param min_read_length numeric, minimum read length in H5 output; Default = 10 (set in generate_stats_figs.R from yaml)
#' @param asite_displacement_length data frame with columns `read_length` and `asite_displacement`
#'  default: read_length = c(28, 29, 30), and asite_displacement = c(15, 15, 15).
#' @param colsum_out logical; if true, return summary column of summed a-site lengths; default: TRUE
#'
#' @return numeric vector if colsum_out = TRUE; matrix with number of rows equivalent to number of rows in asite_displacement_length if colsum_out=FALSE
#'
#' @examples
#' reads_pos_length <- GetGeneDatamatrix(gene = "YAL068C", dataset = "vignette", hd_file = "vignette/output/WTnone/WTnone.h5")
#'  # int [1:41, 1:863] 0 0 0 0 0 0 0 0 0 0 ...
#'
#' CalcAsiteFixed(reads_pos_length, min_read_length = 10, asite_displacement_length = data.frame(read_length = c(28, 29, 30), asite_displacement = c(15, 15, 15)), colsum_out = TRUE)
#'  # num [1:863] 0 0 0 0 0 0 0 0 0 0 ...
#'
#' @export
CalcAsiteFixed <- function(reads_pos_length, min_read_length,
                           asite_displacement_length = data.frame(
                             read_length = c(28L, 29L, 30L),
                             asite_displacement = c(15L, 15L, 15L)
                           ),
                           colsum_out = TRUE) {
  npos <- ncol(reads_pos_length)
  Asite_counts_bylength <-
    purrr::map2(
      asite_displacement_length$read_length, asite_displacement_length$asite_displacement,
      function(read_length, asite_displacement) {
        CalcAsiteFixedOneLength(
          reads_pos_length,
          min_read_length,
          read_length,
          asite_displacement
        )
      }
    )
  if (colsum_out) {
    Asite_counts <- purrr::reduce(Asite_counts_bylength, `+`)
    return(Asite_counts)
  } else {
    # this has only as many columns as asite_displacement_length,
    # probably LESS than data_mat
    Asite_counts_bylengthmat <- unlist(Asite_counts_bylength) %>%
      matrix(ncol = npos, byrow = TRUE)
    return(Asite_counts_bylengthmat)
    # nrow(Asite_counts_bylengthmat) represents the rows in read_length column of asite_displacement_length.
    # TODO: perhaps add row naming to make this ^ clear?
  }
}



createCountDataFrame <- function(gene_name,start,end,codons,h5_file,dataset,asite_displacement_length,min_read_length,frame0_only=F)
{
  count.matrix <- GetGeneDatamatrix(gene_name,dataset,h5_file)
  counts.by.nt <- CalcAsiteFixed(count.matrix,min_read_length,asite_displacement_length)
  if (frame0_only)
  { 
    frame0 <- seq(start,end,3)
    counts.by.codon <- counts.by.nt[frame0]
  } else{
    counts.by.nt <- counts.by.nt[start:end]
    counts.by.codon <- zoo::rollapply(counts.by.nt,3,sum,align="left",by=3)
  }
  if(length(counts.by.codon) != length(codons))
  {
    stop("Error: List of codons in sequence does not match length of RFP counts")
  }
  codon.count.df <- data.frame(gene=rep(gene_name,length(counts.by.codon)),
                               Position=seq(1,length(counts.by.codon)),
                               Codon = codons,
                               RFPCount=counts.by.codon
  )
  codon.count.df <- codon.count.df %>% filter(Position > 1) %>% mutate(Position = Position - 1)
  return(codon.count.df)
}

getCodons <- function(sequence,start,end)
{
  return(as.character(Biostrings::codons(sequence[start:end])))
}

#' readGFFAsDf(): Read GFF file as a tibble (nicer dataframe)
#'
#' Read in positions of all reatures in GFF format and convert to tibble data frame
#'
#' @param orf_gff_file A filepath to a riboviz generated GFF2/GFF3 annotation file
#'
#' @return Tidy data frame (tibble) of GFF data from GFF file
#'
#' @examples
#' readGFFAsDf(orf_gff_file="vignette/input/yeast_YAL_CDS_w_250utrs.gff3")
#'
#' @export
readGFFAsDf <- purrr::compose(
  rtracklayer::readGFFAsGRanges,
  data.frame,
  as_tibble,
  .dir = "forward" # functions called from left to right
)



cds.seq.file <- "~/example-datasets/fungi/saccharomyces/annotation/Saccharomyces_cerevisiae_yeast_CDS_w_250utrs.fa"
asite.displacement.length.file <- "../00_data/02_riboseq_asite/chou_asite_offset.txt"
cds.seq <- readDNAStringSet(cds.seq.file)

gff_df <- readGFFAsDf("~/example-datasets/fungi/saccharomyces/annotation/Saccharomyces_cerevisiae_yeast_CDS_w_250utrs.gff3") %>%
  filter(type == "CDS")
start <- gff_df$start
end <- gff_df$end - 3 # remove stop codon
gene_names <- as.character(unique(gff_df$seqnames))

h5_file <- "/home/copea1/Public_sequencing/Ribo_seq/Fungi/Scerevisiae/Chou_etal_2017_Mol_Cell/output_trimmed/elp1D_2/elp1D_2.h5"
dataset <- "C-Sc_2017"
count_threshold <- 64
min_read_length <- 10
asite_displacement_length <- read_tsv(asite.displacement.length.file,comment="#")

count_df <- purrr::pmap(list(cds.seq,gene_names,start,end),function(gene,gene_name,start,end)
  {
    codons <- getCodons(gene,start,end)
    createCountDataFrame(gene_name = gene_name,
                         start = start,
                         end = end,
                         codons = codons,
                         h5_file = h5_file,
                         dataset = dataset,
                         asite_displacement_length = asite_displacement_length,
                         min_read_length = min_read_length,
                         frame0_only=F)
}) %>% bind_rows()

write_csv(count_df,"../00_data/00_panse_input/00_unfiltered_genes/chou_etal_2019_elp1D_2_all_genes.csv")



