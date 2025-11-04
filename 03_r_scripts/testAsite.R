library(tidyverse)
library(rhdf5)
library(Biostrings)
library(rtracklayer)

source(file.path("/home","copea1","riboviz","rscripts","provenance.R"))
source(file.path("/home","copea1","riboviz","rscripts", "read_count_functions.R"))
source(file.path("/home","copea1","riboviz","rscripts", "stats_figs_block_functions.R"))

CalculateCodonSpecificRibosomeDensity_local <- function(t_rna_file, codon_positions_file, gene_names, hd_file, dataset, gff_df, count_threshold, a_site_displacement){
	
	
	trna <- read_tsv(t_rna_file) 
	load(codon_positions_file) # Position of codons in each gene (numbering ignores first 200 codons)
	# Reads in an object named "codon_pos"
	
	# Pull CDS start and end positions
	gff_df_cds <- gff_df %>% filter(type=="CDS")
	
	start_pos <- gff_df_cds$start
	end_pos <- gff_df_cds$end
	names(start_pos) <- gff_df_cds$seqnames
	names(end_pos) <- gff_df_cds$ seqnames
	
	a_site_displacement_min_read_length <- a_site_displacement %>% filter(read_length >= min_read_length)
	
	out <- lapply(gene_names, function(gene) {
		# From "Position specific distribution of reads" plot
		GetCodonPositionReads(gene, dataset, 
													hd_file = hd_file, 
													left = start_pos[gene], 
													right = end_pos[gene], 
													min_read_length = min_read_length, 
													a_site_displacement = a_site_displacement_min_read_length)
		
	}) # Get codon-based position-specific reads for each gene
	names(out) <- gene_names
	
	gene_len <- sapply(out, length) # Calculate gene length in codons
	out <- out[gene_len > 201] # Ignore genes with <=200 sense codons
	
	trim_out <- lapply(out, function(x) {
		x[201:(length(x) - 1)]
	}) # Trim first 200 codons and stop codon from each gene
	read_counts_trim <- sapply(trim_out, sum) # Calculate read counts in trimmed genes
	trim_out <- trim_out[read_counts_trim >= count_threshold] # Ignore genes with fewer than count_threshold mapped reads
	
	norm_out <- lapply(trim_out, function(x) {
		x / mean(x)
	}) # Normalize reads in each gene by their mean
	
	# Calculate codon-specific mean ribosome-densities at A/P/E sites of the mapped reads
	a_mn <- sapply(names(codon_pos), function(codon) {
		mean(unlist(apply(codon_pos[[codon]], 1, function(a) {
			pos <- as.numeric(a[2])
			norm_out[[a[1]]][pos]
		})), na.rm = T)
	})
	p_mn <- sapply(names(codon_pos), function(codon) {
		mean(unlist(apply(codon_pos[[codon]], 1, function(a) {
			pos <- as.numeric(a[2]) + 1 ## original code has +1
			
			norm_out[[a[1]]][pos]
		})), na.rm = T)
	})
	e_mn <- sapply(names(codon_pos), function(codon) {
		mean(unlist(apply(codon_pos[[codon]], 1, function(a) {
			pos <- as.numeric(a[2]) + 2 ## original code has +2
			
			norm_out[[a[1]]][pos]
		})), na.rm = T)
	})
	
	# Sort the values
	A <- a_mn[order(names(codon_pos))]
	P <- p_mn[order(names(codon_pos))]
	E <- e_mn[order(names(codon_pos))]
	trna <- trna[order(trna$Codon),]
	cod_dens_tRNA_data <- cbind(trna, A, P, E)
	return(cod_dens_tRNA_data)
	
} # end of CalculateCodonSpecificRibosomeDensity() definitionlateCodonSpecificRibosomeDensity() definition

readGFFAsDf <- purrr::compose(
	rtracklayer::readGFFAsGRanges,
	data.frame,
	as_tibble,
	.dir = "forward" # functions called from left to right
)



cds.seq.file <- "~/example-datasets/fungi/saccharomyces/annotation/Saccharomyces_cerevisiae_yeast_CDS_w_250utrs.fa"
asite.displacement.length.file <- "../00_data/02_riboseq_asite/wu_asite_offset.txt"
t_rna_file <- "~/riboviz/data/yeast_tRNAs.tsv"
codon_positions_file <- "~/riboviz/data/yeast_codon_pos_i200.RData"


cds.seq <- readDNAStringSet(cds.seq.file)
gff_df <- readGFFAsDf("~/example-datasets/fungi/saccharomyces/annotation/Saccharomyces_cerevisiae_yeast_CDS_w_250utrs.gff3") %>%
	filter(type == "CDS")
start <- gff_df$start
end <- gff_df$end - 3 # remove stop codon
gene_names <- as.character(unique(gff_df$seqnames))

hd_file <- "/nobackup/rokaslab/copea1/Public_sequencing/Ribo_seq/Fungi/Scerevisiae/Wu_etal_2019_Mol_Cell/output/CHX_1/CHX_1.h5"
dataset <- "Wu_etal_2019"
count_threshold <- 64
min_read_length <- 10
a_site_displacement <- read_tsv(asite.displacement.length.file,comment="#")


new.wt<-CalculateCodonSpecificRibosomeDensity_local(t_rna_file, 
																			codon_positions_file, 
																			gene_names, 
																			hd_file, 
																			dataset, 
																			gff_df, 
																			count_threshold, 
																			a_site_displacement)
	
	