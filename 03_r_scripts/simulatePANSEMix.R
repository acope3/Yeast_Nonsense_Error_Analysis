library(argparse)

rm(list=ls())


parser <- ArgumentParser()
parser$add_argument("-i","--input",type="character",default="./")
parser$add_argument("-o","--output",type="character",default="./")
parser$add_argument("--alpha",type="character",default="pa_Alpha.csv")
parser$add_argument("--lambda",type="character",default="pa_LambdaPrime.csv")
parser$add_argument("--nserate",type="character",default="nse_rate_10_5.csv")
parser$add_argument("--phi",type="character",default="Weinberg_2016/scer_roc_phi_locus_for_weinberg_4000_genes_500_counts_min.csv")
parser$add_argument("--footprint_count",type="integer",default=NULL)
parser$add_argument("--resample_phi",action="store_true")
parser$add_argument("--development",help="Use this flag to run a developmental version of AnaCoDa",type="character",default=NULL)



args <- parser$parse_args()
alpha <- args$alpha
input <- args$input
output <- args$output
lambda <- args$lambda
nserate <- args$nserate
phi <- args$phi
footprint_count <- args$footprint_count
resample_phi <- args$resample_phi
dev <- args$development


if (!is.null(dev))
{
  library(AnaCoDa,lib.loc=dev)
} else {
  library(AnaCoDa,lib.loc = "~/AnaCoDa_installs/PANSE_allow_different_wt_sim/")
}


dir.create(output)
cmd <- paste(commandArgs(T),collapse = " ")
cmd <- paste("Rscript --vanilla R_scripts/simulatePANSEMix.R",cmd,sep=" ")
readme <- paste("Simulation was performed with the command\n```\n",cmd,"\n```",sep="")
write(readme,file.path(output,"README.md"))


fasta.files <- input
numMixtures <- 1
genome <- initializeGenomeObject(file=fasta.files,fasta=F,positional=T)
#initialize parameter object

rfp <- read.table(fasta.files,sep=",",header=T,stringsAsFactors = F)

numElongationMixtures <- length(unique(rfp$Mixture))

sphi_init <- rep(1,numMixtures)
with.phi <- F
mixDef <- "allUnique"
size <- length(genome)
index <- c(1:size)

geneAssignment <- c(rep(1, length(genome))) 

segment_exp <- read.table(file=phi,sep=",",header=TRUE)
##segment_exp[,phi_type] <- segment_exp[,phi_type]/mean(segment_exp[,phi_type])
init_phi <-segment_exp[,2]


sphi_init <- sd(log(init_phi))

alpha.files <- unlist(strsplit(alpha,split=","))
lambda.files <- unlist(strsplit(lambda,split=","))
nse.files <- unlist(strsplit(nserate,split=","))

if (resample_phi)
{
	mean.phi <- 0
	while(abs(mean.phi - 1) > 0.05)
	{
		rand_phi <- rlnorm(length(init_phi),meanlog=-(sphi_init^2)/2,sdlog=sphi_init)
		mean.phi <- mean(rand_phi)
	}
	rand_phi <- sort(rand_phi)
	init_phi <- rand_phi[rank(init_phi)]
}

if (is.null(footprint_count))
{
	footprint_count <- sum(rfp[,"RFPCount"])
}
genome$setSumRFP(footprint_count)

parameter <- initializeParameterObject(genome,model="PANSE",sphi_init,numMixtures, geneAssignment, split.serine = FALSE, mixture.definition = mixDef, initial.expression.values = init_phi,numElongationMixtures = numElongationMixtures)
parameter$initMutationSelectionCategories(c(alpha.files), length(alpha.files), "Alpha")
parameter$initMutationSelectionCategories(c(lambda.files), length(lambda.files),"LambdaPrime")
parameter$initMutationSelectionCategories(c(nse.files), length(nse.files),"NSERate")
model <- initializeModelObject(parameter, "PANSE",rfp.count.column=1)
model$simulateGenome(genome)


output.file <- paste(output,"simulated_rfp.csv",sep="/")

genome$writeRFPData(output.file,T)


sim.rfp <- read.table(output.file,sep=",",header=T,stringsAsFactors = F)
sim.rfp <- sim.rfp[which(!sim.rfp$Codon %in% c("TAA","TAG","TGA")),]
write.table(sim.rfp,output.file,sep=",",row.names=F,col.names=T,quote=F)
segment_exp[,2] <- init_phi
write.table(segment_exp[,c(1,2)],paste(output,"phi.csv",sep="/"),sep=",",row.names=F,col.names=T,quote=F)

