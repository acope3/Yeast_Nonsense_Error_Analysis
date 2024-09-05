#!/usr/bin/env Rscript

library(profmem)
library(argparse)
rm(list=ls())

parser <- ArgumentParser()
parser$add_argument("-i","--input",help="File storing RFP data",type="character",default="./")
parser$add_argument("-o","--output",help="Directory of where to put results. Will automatically generate lowest-level directory if not already generated.",type="character",default="./")
parser$add_argument("-d","--div",help="Number of steps to diverge from starting values. Will be applied at beginning of each run, with the exception of the last.",type="integer",default=0)
parser$add_argument("-s","--samp",help="Number of samples",type="integer",default=5000)
parser$add_argument("-a","--adapt",help="Adaptive Width",type="integer",default=50)
parser$add_argument("-t","--thin",help="Thinning value",type="integer",default=20)
parser$add_argument("-n","--threads",help="Number of threads to use for MCMC",type="integer",default=1)
parser$add_argument("--alpha",help="Initial Alpha values. Assumes csv format with Codons in first column and Alpha values in second column",type="character",default=NULL)
parser$add_argument("--lambda",help="Initial Lambda Prime values. Assumes csv format with Codons in first column and Lambda Prime values in second column",type="character",default=NULL)
parser$add_argument("--nserate",help="Initial NSE Rate values. Assumes csv format with Codons in first column and NSE Rate values in second column",type="character",default=NULL)
parser$add_argument("--phi",help="Initial Phi values. Assumes csv format with Gene IDs in first column and Phi values in second column", type="character")
parser$add_argument("--normalize_phi",help="Normalize phi such that the initial values have arithmetic mean of 1 on the natural scale.",action="store_true")
parser$add_argument("--est_csp",help="Use this flag to indicate estimation of CSP. Otherwise, CSP will not be estimated",action="store_true")
parser$add_argument("--est_phi",help="Use this flag to indicate estimation of Phi. Otherwise, Phi will not be estimated",action="store_true")
parser$add_argument("--est_hyp",help="Use this flag to indicate estimation of Hyperparameters. Otherwise, Hyperparameters will not be estimated",action="store_true")
parser$add_argument("--max_num_runs",help="Max number of runs to do. Note that each fitting will consists of at least two runs. The first can serve as a burn-in.",type="integer",default = 6)
parser$add_argument("--fix_nse",help="Use this flag to fix NSERate at starting value.",action="store_true")
parser$add_argument("--fix_alpha",help="Use this flag to fix alpha at starting value",action="store_true")
parser$add_argument("--fix_lambda",help="Use this flag to fix lambda^prime at starting value.",action="store_true")
parser$add_argument("--fix_sphi",help="Use this flag to fix s_phi at starting value.",action="store_true")
parser$add_argument("--init_z",type="double",default=NULL)
parser$add_argument("--share_nse",action="store_true")
parser$add_argument("--prior_type",type="character",default="Natural_Uniform")
parser$add_argument("--nserate_uniform_lower_limit",type="double",default=1e-100)
parser$add_argument("--nserate_uniform_upper_limit",type="double",default=1e-1)
parser$add_argument("--nserate_exponential_mean",type="double",default=25000)
parser$add_argument("--mixture_definition",type="character",help="Must be allUnique (default), elongationShared, or nseShared.",default="allUnique")
parser$add_argument("--fix_z",help="Use this flag to fix s_phi at starting value.",action="store_true")
parser$add_argument("--phi_column",type="integer",default=2)
parser$add_argument("--with_phi",action="store_true")
parser$add_argument("--observed_phi",type="character")
parser$add_argument("--s_epsilon",type="double",default=0.05)
parser$add_argument("--s_phi",type="double",default=NULL)
parser$add_argument("--dataset",type="character",default="RFP")
parser$add_argument("--restart_file",type="character",default=NULL)
parser$add_argument("--development",help="Use this flag to run a developmental version of AnaCoDa",type="character",default=NULL)


args <- parser$parse_args()
div <- args$div
input <- args$input
directory <- args$output
thinning <- args$thin
adaptiveWidth <- args$adapt
samples <- args$samp
num_threads <- args$threads
phi.files <- args$phi
normalize.phi <- args$normalize_phi
alpha <- args$alpha
lambda <- args$lambda
nserate <- args$nserate
mixture.definition <- args$mixture_definition
est.phi <- args$est_phi
est.csp <- args$est_csp
est.hyp <- args$est_hyp
max_num_runs <- args$max_num_runs
fix_nse <- args$fix_nse
fix_alpha <- args$fix_alpha
fix_lambda <- args$fix_lambda
fix_sphi <- args$fix_sphi
share_nse <- args$share_nse
prior.type <- args$prior_type
nserate.uniform.lower.limit <- args$nserate_uniform_lower_limit
nserate.uniform.upper.limit <- args$nserate_uniform_upper_limit
nserate.exponential.mean <- args$nserate_exponential_mean
init_z <- args$init_z
fix_z <- args$fix_z
with.phi <- args$with_phi
observed.phi <- args$observed_phi
restart.file <- args$restart_file
s.epsilon <- args$s_epsilon
s.phi <- args$s_phi
dataset <- args$dataset
phi_column <- args$phi_column
dev <- args$development

print(args)

if (!is.null(dev))
{
  library(AnaCoDa,lib.loc=dev)
} else {
  library(AnaCoDa)
}
createParameterOutput <- function(parameter,numMixtures,samples,mixture.labels,samples.percent.keep=1,relative.to.optimal.codon=F,report.original.ref=T)
{
  for (i in 1:numMixtures)
  {
    getCSPEstimates(parameter,paste(dir_name,"Parameter_est",mixture.labels[i],sep="/"),i,samples*samples.percent.keep,relative.to.optimal.codon=relative.to.optimal.codon,report.original.ref = report.original.ref)
    getCSPEstimates(parameter,paste(dir_name,"Parameter_est",mixture.labels[i],sep="/"),i,samples*samples.percent.keep,relative.to.optimal.codon=relative.to.optimal.codon,report.original.ref = report.original.ref,log.scale=T)  
  }
}

createTracePlots <- function(trace, model,genome,numMixtures,samples,mixture.labels,samples.percent.keep=1)
{
  for (i in 1:numMixtures)
  {
    plot(trace, what = "Alpha", mixture = i)
    plot(trace, what="Lambda",mixture =i)
    plot(trace, what = "MeanWaitingTime", mixture = i)
    plot(trace, what = "VarWaitingTime", mixture = i)
    plot(trace,what="NSERate", mixture = i)
    plot(trace,what="NSERate",log.10.scale=T, mixture = i)
    plot(trace,what="NSEProb", mixture = i)
    plot(trace,what="NSEProb",log.10.scale=T, mixture = i)
    
  }
  plot(trace,what="AcceptanceRatio")
}

createAcceptanceRateTrace <- function(csp.accept.trace,param="Elongation",percent.to.keep=0.5)
{
  done.adapt <- TRUE
  codon.list <- codons()
  for(i in 1:61)
  {
    accept.trace <- csp.accept.trace[[i]]
    len <- length(accept.trace)
    mean.acceptance <- mean(accept.trace[(len-len*percent.to.keep):len])
    if (mean.acceptance < 0.1 || mean.acceptance > 0.44)
    { 
      done.adapt <- FALSE
    }
    plot(accept.trace,main=paste0(param," Acceptace Rate for ",codon.list[i]),xlab="Samples",ylab="Acceptance Rate",type="l")
  }
  return(done.adapt)
}


estimateStartingZ <- function(genome,alpha,lp,phi)
{
  alpha.files <- unlist(strsplit(alpha,split=","))
  lambda.files <- unlist(strsplit(lp,split=","))
  alpha.val <- read.table(alpha.files[1],sep=",",header=T,stringsAsFactors=F)
  lp.val <- read.table(lambda.files[1],sep=",",header=T,stringsAsFactors=F)
  wait.times <- alpha.val[,2]/lp.val[,2]
  codon.counts.per.gene <- getCodonCounts(genome)
  codon.counts.per.gene <- data.matrix(codon.counts.per.gene[2:62]) #exclude gene name and stop codons
  elong.time.mrna <- codon.counts.per.gene %*% wait.times
  total <- elong.time.mrna * phi
  z <- sum(total)
  return(z)
}
dir.create(directory)

cmd <- paste(commandArgs(T),collapse = " ")
cmd <- paste("Rscript --vanilla R_scripts/runPANSEMixModel.R",cmd,sep=" ")
readme <- paste("Model was run with bash command\n```\n",cmd,"\n```",sep="")
write(readme,file.path(directory,"README.md"))

fasta.files <- input
numMixtures <- 1

## Note: writing a for loop to deal with all mixtures (1 - n.mixtures) is tricky.
## Part of the issue is the appending of the object defined in the command and the assignment of the outputT


if (with.phi)
{
  genome <- initializeGenomeObject(file=fasta.files,fasta=F,match.expression.by.id=F,observed.expression.file=observed.phi)
} else{
  genome <- initializeGenomeObject(file=fasta.files,fasta=F)

}


print(length(genome))
cat("Genome loaded\n")
#initialize parameter object



mixDef <- mixture.definition
percent.to.keep <- 0.5
size <- length(genome)
cat(size,"\n")
index <- c(1:size)

geneAssignment <- c(rep(1, length(genome))) 

if (!is.null(phi.files))
{
  segment_exp <- read.table(file=phi.files,sep=",",header=TRUE)
  init_phi <-segment_exp[,phi_column]
  if (normalize.phi)
  {
    init_phi <- init_phi/mean(init_phi)
  }
  obs_sphi_init <- rep(sd(log(init_phi)),numMixtures)
  if(length(genome) != length(init_phi)){
  stop("length(genomeObj) != length(init_phi), but it should.")
  }else{
    print("Initial Phi values successfully files loaded:");
  }

} else{
  init_phi <-NULL
}

if (!is.null(s.phi))
{
  sphi_init <- rep(s.phi,numMixtures)
} else if (!is.null(phi.files)){
  sphi_init <- rep(obs_sphi_init,numMixtures)
}else{
  sphi_init <- rep(1.5,numMixtures)
}


rfp <- read.table(fasta.files,sep=",",header=T,stringsAsFactors = F)
numElongationMixtures <- length(unique(rfp$Mixture^2)) ## will handle negatives
mixture.labels <- paste(seq(numElongationMixtures),dataset,sep="_")
                        
if (is.null(init_z))
{
  if (!is.null(alpha) && !is.null(lambda) && !is.null(init_phi))
  {
    init_z <- estimateStartingZ(genome,alpha,lambda,init_phi)
  } else {
    init_z <- 100000
  }
}
meta.title <- unlist(strsplit(directory,"/"))
meta.title <- meta.title[length(meta.title)]


done <- FALSE
run_number <- 1
while(run_number <= max_num_runs)
{
  
  if (run_number == 1)
  {
    if (is.null(restart.file))
    {
      parameter <- initializeParameterObject(genome,model="PANSE",sphi_init,numMixtures, geneAssignment, split.serine = FALSE, mixture.definition = mixDef, initial.expression.values = init_phi,init.partition.function=init_z, init.sepsilon = s.epsilon, numElongationMixtures = numElongationMixtures)
        
      if (!is.null(alpha))
      {
        alpha.files <- unlist(strsplit(alpha,split=","))
        parameter$initMutationSelectionCategories(alpha.files, length(alpha.files), "Alpha")
      }
      if (!is.null(lambda))
      {
        lambda.files <- unlist(strsplit(lambda,split=","))
        parameter$initMutationSelectionCategories(lambda.files, length(lambda.files), "LambdaPrime")
      }
      if (!is.null(nserate))
      {
        nse.files <- unlist(strsplit(nserate,split=","))
        parameter$initMutationSelectionCategories(nse.files, length(nse.files), "NSERate")
      }
      dir_name <- paste0(directory,"/restart_",run_number)
      steps.to.adapt <- (samples*thinning)*percent.to.keep
      if (run_number == max_num_runs)
      {
        dir_name <- paste0(directory,"/final_restart")
      }
    } else{
      print("Initialize from Restart File")
      previous <- stringr::str_extract(pattern="restart_[0-9]+",string=restart.file)
      if (!is.na(previous))
      {
        run_number <- as.numeric(stringr::str_extract(pattern="[0-9]+",string=previous)) + 1
      }
      parameter<-initializeParameterObject(init.with.restart.file = restart.file,model="PANSE")
      steps.to.adapt <- (samples*thinning)*percent.to.keep
      dir_name <- paste0(directory,"/restart_",run_number)
      
      if (run_number == max_num_runs)
      {
        steps.to.adapt <- 0
        dir_name <- paste0(directory,"/final_restart")
      }
    }
    
  } else if (run_number == max_num_runs){
    steps.to.adapt <- 0
    parameter<-initializeParameterObject(init.with.restart.file = paste(dir_name,"Restart_files/rstartFile.rst_final",sep="/"),model="PANSE")
    dir_name <- paste0(directory,"/final_restart")
  } else{

    parameter<-initializeParameterObject(init.with.restart.file = paste(dir_name,"Restart_files/rstartFile.rst_final",sep="/"),model="PANSE")
    dir_name <- paste0(directory,"/restart_",run_number)
    steps.to.adapt <- (samples*thinning)*percent.to.keep
    
  }
  if (fix_nse)
  {
    parameter$fixNSERate()
  }
  if (fix_alpha)
  {
    parameter$fixAlpha()
  }
  if (fix_lambda)
  {
    parameter$fixLambdaPrime()
  }
  if(fix_sphi)
  {
    parameter$fixSphi()
  }
  if(fix_z)
  {
    parameter$fixZ()
  }
  if (share_nse)
  {
    parameter$shareNSERate()
  }
  
  print(paste0("Creating",dir_name))
  dir.create(dir_name)
  dir.create(paste(dir_name,"Graphs",sep="/"))
  dir.create(paste(dir_name,"Restart_files",sep="/"))
  dir.create(paste(dir_name,"Parameter_est",sep="/"))
  dir.create(paste(dir_name,"R_objects",sep="/"))
  
  mcmc <- initializeMCMCObject(samples=samples, thinning=thinning, adaptive.width=adaptiveWidth,
                                 est.expression=est.phi, est.csp=est.csp, est.hyper=est.hyp,est.mix=FALSE)
  
  mcmc$setStepsToAdapt(steps.to.adapt)
  
  
  model <- initializeModelObject(parameter, "PANSE", with.phi,fix.observation.noise=T)
  model$setNSERatePriorDistribution(prior.type,
                                    nserate.uniform.lower.limit,
                                    nserate.uniform.upper.limit,
                                    nserate.exponential.mean)
  
  setRestartSettings(mcmc, paste(dir_name,"Restart_files/rstartFile.rst",sep="/"), adaptiveWidth, F)
  sys.runtime <- system.time(
    runMCMC(mcmc, genome, model, num_threads,div=div)
  )
  sys.runtime <- data.frame(Value=names(sys.runtime),Time=as.vector(sys.runtime))
  write.table(sys.runtime,file=paste(dir_name,"mcmc_runtime.csv",sep="/"),sep=",",col.names = T,row.names = T,quote=F)
  
  writeParameterObject(parameter,paste(dir_name,"R_objects/parameter.Rda",sep="/"))
  writeMCMCObject(mcmc,file=paste(dir_name,"R_objects/mcmc.Rda",sep="/"))
  
  createParameterOutput(parameter = parameter,numMixtures = numElongationMixtures,mixture.labels = mixture.labels,samples = samples,samples.percent.keep = percent.to.keep,relative.to.optimal.codon = F,report.original.ref = T)
  expressionValues <- getExpressionEstimates(parameter,c(1:size),samples*percent.to.keep,genome=genome)
  write.table(expressionValues,file=paste(dir_name,"Parameter_est/gene_expression.txt",sep="/"),sep=",",col.names = T,quote = F,row.names = F)
  
  
  
  trace <- parameter$getTraceObject()
  if (est.csp)
  {
    pdf(paste(dir_name,"Graphs/CSP_traces.pdf",sep="/"), width = 11, height = 12,title=paste0("CSP_traces_",meta.title,"_restart_",run_number,".pdf"))
    createTracePlots(trace=trace,model=model,genome=genome,numMixtures=numElongationMixtures,samples=samples,samples.percent.keep = 1,mixture.labels = mixture.labels)
    nse.accept.trace <- trace$getNseRateSpecificAcceptanceRateTrace()
    nse.done.adapt <- createAcceptanceRateTrace(nse.accept.trace,"NSE Rate",percent.to.keep) 
    dev.off()
  }
  #plots different aspects of trace
  
  param.conv <- TRUE
  pdf(paste(dir_name,"Graphs/mcmc_traces.pdf",sep="/"),title=paste0("mcmc_traces_",meta.title,"_restart_",run_number,".pdf"))
  plot(mcmc,what = "LogPosterior")
  plot(mcmc,what = "LogLikelihood")
  plot(trace, what = "ExpectedPhi")
  if (est.hyp)
  {
    if (!fix_sphi)
    {
      plot(trace,what="Sphi")
      plot(trace,what="Mphi")
    }
    if (!fix_z)
    {
      plot(trace,what="PartitionFunction")
    }
    if (with.phi)
    {
      plot(trace,what="Sepsilon")
    }
  }
  if (est.csp)
  {
    if (!fix_alpha)
    {
      acfCSP(parameter,csp="Alpha",numMixtures = numMixtures,samples=samples*percent.to.keep)
    }
    if (!fix_lambda)
    {
      acfCSP(parameter,csp="LambdaPrime",numMixtures = numMixtures,samples=samples*percent.to.keep)
    }
    if (!fix_nse)
    {
      acfCSP(parameter,csp="NSERate",numMixtures = numMixtures,samples=samples*percent.to.keep)
    }
  }
  dev.off()
  
  if (est.csp)
  {
    for (i in 1:numMixtures)
    {
      param.diag<-convergence.test(trace,samples=samples*percent.to.keep,thin = thinning,what="Alpha",mixture=i,frac1=0.1)
      z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
      write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_alpha_",i,".txt"),ncolumns = 1)
    }
    
    
    for (i in 1:numMixtures)
    {
      param.diag<-convergence.test(trace,samples=samples*percent.to.keep,thin = thinning,what="LambdaPrime",mixture=i,frac1=0.1)
      z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
      write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_lambda_",i,".txt"),ncolumns = 1)
    }
    
    for (i in 1:numMixtures)
    {
      param.diag<-convergence.test(trace,samples=samples*percent.to.keep,thin = thinning,what="NSERate",mixture=i,frac1=0.1)
      z.scores <- param.diag$z[which(abs(param.diag$z) > 1.96)]
      write(param.diag$z,paste0(dir_name,"/Parameter_est/convergence_nserate_",i,".txt"),ncolumns = 1)
    }
  }
  rm(parameter)
  rm(trace)
  rm(model)
  run_number <- run_number + 1
}

