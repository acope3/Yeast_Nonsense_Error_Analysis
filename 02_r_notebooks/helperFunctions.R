library(tidyverse)
library(ggpubr)
library(AnaCoDa)
library(cowplot)
library(ggrepel)
library(bayestestR)

createPattern <- function(number,string)
{
  return(paste0(number,"_.*",string))
}

getParameterDataFrames <- function(filepath,category=1)
{
  alpha.files <- list.files(file.path(filepath,"Parameter_est"),pattern = createPattern(category,"_Alpha.csv"),full.names = T)
  lambda.files <- list.files(file.path(filepath,"Parameter_est"),pattern = createPattern(category,"_Lambda(_Prime){0,1}.csv"),full.names = T)
  nserate.files <- list.files(file.path(filepath,"Parameter_est"),pattern = createPattern(category,"_NSERate.csv"),full.names = T)
  log.nserate.files <- list.files(file.path(filepath,"Parameter_est"),pattern = createPattern(category,"_NSERate_log_scale.csv"),full.names = T)
  nseprob.files <- list.files(file.path(filepath,"Parameter_est"),pattern = createPattern(category,"_NSEProb.csv"),full.names = T)
  log.nseprob.files <- list.files(file.path(filepath,"Parameter_est"),pattern = createPattern(category,"_NSEProb_log_scale.csv"),full.names = T)
  wait.time.files <- list.files(file.path(filepath,"Parameter_est"),pattern = createPattern(category,"_WaitingTime.csv"),full.names = T)
  phi.files <- list.files(file.path(filepath,"Parameter_est"),pattern = "gene_expression.txt",full.names = T)
  
  
  alpha <- read_csv(alpha.files,col_types = cols())
  lambda <- read_csv(lambda.files,col_types = cols())
  if (length(wait.time.files) == 0)
  {
    wait.time <- data.frame(AA=alpha$AA,Codon=alpha$Codon,Mean=alpha$Mean/lambda$Mean)
  } else{
    wait.time <- read_csv(wait.time.files,col_types = cols())
  }
  nserate <- read_csv(nserate.files,col_types = cols())
  log.nserate <- read_csv(log.nserate.files,col_types = cols())
  nseprob <- read_csv(nseprob.files,col_types = cols())
  log.nseprob <- read_csv(log.nseprob.files,col_types = cols())
  phi <- read_csv(phi.files,col_types = cols())
  
  params <- list(Alpha = alpha,
                 Lambda = lambda,
                 Wait.Time = wait.time,
                 NSERate = nserate,
                 Log.NSERate = log.nserate,
                 NSEProb = nseprob,
                 Log.NSEProb = log.nseprob,
                 Phi = phi)
  return(params)
}

nucleotideDiff <- function(codon1,codon2)
{
  codon1 <- unlist(str_split(codon1,""))
  codon2 <- unlist(str_split(codon2,""))
  diff.nuc.pos <- codon1 != codon2
  diff.nuc <- paste0(sort(c(codon1[diff.nuc.pos],codon2[diff.nuc.pos])),collapse = "")
  return(diff.nuc)
}


stopCodonNeighbor <- function(codon,by.pos=0,by.transition = F,by.transversion=F)
{
  if (by.transition == T && by.transversion == T)
  {
    stop("Only by.transition or by.transversion should be true, or neither should be true")
  }
  nucleotides <- unlist(str_split(codon,""))
  first.pos <- str_c(c("A","T","C","G"),nucleotides[2],nucleotides[3])
  second.pos <- str_c(nucleotides[1],c("A","T","C","G"),nucleotides[3])
  third.pos <- str_c(nucleotides[1],nucleotides[2],c("A","T","C","G"))
  if (by.pos == 1)
  {
    neighbors <- first.pos
  } else if (by.pos == 2)
  {
    neighbors <- second.pos
  } else if (by.pos == 3)
  {
    neighbors <- third.pos
  } else
  {
    neighbors <- c(first.pos,second.pos,third.pos)
  }
  num.neighbors <- 0
  if ("TAG" %in% neighbors) 
  {
    diff.nuc <- nucleotideDiff(codon,"TAG")
    if (by.transition && diff.nuc %in% c("AG","CT"))
    {
      num.neighbors <- num.neighbors + 1
    } else if (by.transversion && diff.nuc %in% c("AC","AT","CG","GT")) {
      num.neighbors <- num.neighbors + 1
    } else if (!by.transition && !by.transversion){
      num.neighbors <- num.neighbors + 1
    }
  } 
  if ("TAA" %in% neighbors) 
  {
    diff.nuc <- nucleotideDiff(codon,"TAA")
    if (by.transition && diff.nuc %in% c("AG","CT"))
    {
      num.neighbors <- num.neighbors + 1
    } else if (by.transversion && diff.nuc %in% c("AC","AT","CG","GT")) {
      num.neighbors <- num.neighbors + 1
    } else if (!by.transition && !by.transversion){
      num.neighbors <- num.neighbors + 1
    }
  } 
  if ("TGA" %in% neighbors) 
  {
    diff.nuc <- nucleotideDiff(codon,"TGA")
    if (by.transition && diff.nuc %in% c("AG","CT"))
    {
      num.neighbors <- num.neighbors + 1
    } else if (by.transversion && diff.nuc %in% c("AC","AT","CG","GT")) {
      num.neighbors <- num.neighbors + 1
    } else if (!by.transition && !by.transversion) {
      num.neighbors <- num.neighbors + 1
    } 
  }
  return(as.character(num.neighbors))
}


compareEstimates <- function(df.1,df.2,variable,xlab,ylab,title,color = "AA",log.scale.phi=F)
{
  df.1 <- df.1[[variable]]
  df.2 <- df.2[[variable]]
  if (variable != "Phi")
  {
    merge.df <- df.1 %>% 
      left_join(df.2,by=c("AA","Codon"))
    
    
    merge.df <- merge.df %>% 
      mutate(Stop.Neighbor = as.character(purrr::map_chr(Codon,stopCodonNeighbor)))
    
    p <- ggplot(merge.df,aes(x=Mean.x,y=Mean.y)) 
    if (!is.na(color) && variable != "Phi")
    {
      p <- p + 
        geom_point(aes(color = !!sym(color))) +
        geom_errorbar(aes(ymin=`2.5%.y`,ymax=`97.5%.y`,color = !!sym(color)),alpha=0.2) + 
        geom_errorbarh(aes(xmin=`2.5%.x`,xmax=`97.5%.x`,color = !!sym(color)),alpha=0.2) 
    } else if (is.na(color) && variable != "Phi") {
      p <- p +
        geom_point() +
        geom_errorbar(aes(ymin=`2.5%.y`,ymax=`97.5%.y`),alpha=0.2) + 
        geom_errorbarh(aes(xmin=`2.5%.x`,xmax=`97.5%.x`),alpha=0.2) 
    }
    p <- p +
      geom_abline(intercept=0,slope=1,linetype="dashed") +
      theme_cowplot() +
      xlab(xlab) +
      ylab(ylab) +
      stat_cor(method="spearman",label.sep="\n") +
      #theme(aspect.ratio=1) +
      ggtitle(title)
  } else if (variable == "Phi") {
    merge.df <- df.1 %>% 
      inner_join(df.2,by=c("GeneID"))
    
  
    p <- ggplot(merge.df,aes(x=Mean.x,y=Mean.y)) +
      geom_point() +
      geom_abline(intercept=0,slope=1,linetype="dashed") +
      theme_cowplot() +
      xlab(xlab) +
      ylab(ylab) +
      stat_cor(method="spearman",label.sep="\n") +
      #theme(aspect.ratio=1) +
      ggtitle(title)
    if (log.scale.phi)
    {
      p <- p + scale_x_log10() + scale_y_log10()
    }
  }
  return(p)
}

## `codons` argument puts vertical lines where the codon occurs

plotRFPByPosition<- function(df,gene_id,codons=NULL,real.vs.sim=F,title=NULL)
{
  if (is.null(title)) 
  {
    title <- gene_id
  }
  ## try to suppress dplyr message: Adding missing grouping variables: `GeneID` using `invisible()`
  invisible(real.sim.target <- df %>% 
              filter(GeneID == gene_id)
  )
  real.sim.target.long <- real.sim.target %>%
    pivot_longer(starts_with("RFPCount"),
                 names_to="Data",
                 values_to="RFPCount")
  if (!is.null(codons))
  {
    real.sim.target.long <- real.sim.target.long %>%
      mutate(Codon.of.Interest = ifelse(Codon %in% codons,Position,NA))
  }
  # p <- ggplot(real.sim.target.long,
  #             aes(x=Position,y=RFPCount,color=Data)) +
  #   geom_line(size = 1)
  # if (!is.null(codons))
  # {
  #   positions.codons <- real.sim.target.long %>%
  #     filter(Codon %in% codons) %>%
  #    dplyr::select(Codon, Position) 
  #   p <- p +
  #      geom_vline(data = positions.codons, aes(xintercept=Position))
  # }
  if (real.vs.sim)
  {
    p <- ggplot(real.sim.target.long,aes(x=Position,y=RFPCount,color=Data))
  } else
  {
    p <- ggplot(real.sim.target.long,aes(x=Position,y=RFPCount))
  }
  p <- p +
    geom_line() +
    ggtitle(title) +
    theme_cowplot()
  if (!is.null(codons))
  {
    
    if (length(which(is.na(real.sim.target.long$Codon.of.Interest))) != nrow(real.sim.target.long))
    {
      p <- p + geom_vline(aes(xintercept=Codon.of.Interest),linetype="dashed")
    }
  }
  p <- p + ylab("RFP Count")
  return(p)
}

compareParameterVsPosition <- function(rfp.data,model.fit,
                                       parameter.name="Log.NSEProb",
                                       xlab="",
                                       title="",
                                       include.stop.codon.neighbor=T,
                                       flip.coord = F,
                                       correlation.label.pos = c(-3,0.5))
{
  log.nseprob <- model.fit[[parameter.name]]
  
  rfp.data<- rfp.data %>% 
    group_by(GeneID) %>% 
    mutate(Normalized.Position = Position/n())
  
  codon.avg.pos <- rfp.data %>% 
    ungroup() %>% 
    group_by(Codon) %>% 
    summarize(Total = n(),
              Avg.Pos = mean(Normalized.Position),
              SD.Pos = sd(Normalized.Position),
              Upper = Avg.Pos + 2*SD.Pos/sqrt(Total), 
              Lower = Avg.Pos - 2*SD.Pos/sqrt(Total))# %>%
    #filter(!Codon %in% c("TGG","ATG")) # drop 1-codon AA
  
  codon.avg.pos <- codon.avg.pos %>% 
    left_join(log.nseprob,by="Codon") %>%
    ungroup() %>%
    rowwise() %>%
    mutate(Stop.neighbors=as.character(stopCodonNeighbor(Codon)))
  
  pos.vs.nseprob <- ggplot(codon.avg.pos,aes(x=Mean,y=Avg.Pos,label=Codon)) 
  
  if (include.stop.codon.neighbor){
    pos.vs.nseprob <- pos.vs.nseprob +
      geom_point(aes(color=Stop.neighbors))
  } else {
    pos.vs.nseprob <- pos.vs.nseprob +
      geom_point()
  }
  pos.vs.nseprob <- pos.vs.nseprob +
    geom_errorbar(aes(ymin=Lower,ymax=Upper),alpha=0.25) +
    geom_errorbarh(aes(xmin=`2.5%`,xmax=`97.5%`),alpha=0.25) +
    geom_hline(yintercept = 0.5,linetype="dashed") +
    ylim(c(0.45,0.55)) +
    theme_cowplot() +
    xlab(xlab) +
    ylab("Average Normalized Position in Gene") +
    ggtitle(title)
  if (flip.coord)
  {
    pos.vs.nseprob <- pos.vs.nseprob +
      coord_flip()
  }
  
  pos.vs.nseprob <- pos.vs.nseprob +
    stat_cor(method="spearman",label.sep="\n",
             label.x = correlation.label.pos[1],
             label.y=correlation.label.pos[2]) 
  pos.vs.nseprob
}


calculateMissingNext <- function(rfp)
{
  rfp <- rfp %>% 
    group_by(GeneID) %>% 
    mutate(Next = lead(RFPCount),
           Missing = ifelse(Next == 0, 1, 0))
  missing.next <- rfp %>% 
    group_by(Codon) %>% 
    summarize(Frequency = sum(Missing,na.rm=T)/n())
  return(missing.next)
}

calculateSuccessProbability <- function(df)
{
  df <- df %>% 
    mutate(Prob.Successful.Per.Pos = (1/Wait.Time)/(NSERate + 1/Wait.Time),
           Pr.NSE = NSERate/(NSERate + 1/Wait.Time)) %>%
    group_by(GeneID) %>%
    mutate(Sigma.i = cumprod(Prob.Successful.Per.Pos)) %>%
    ungroup()
  return(df)
}


getPPIWidth <- function(nse.df,lower,upper)
{
  nse.df <- nse.df %>%
    mutate(PPI.Width = !!as.name(upper) - !!as.name(lower))
  return(nse.df)
}

getNSEForReadDepth <- function(model.fit.dir,param.type="Log.NSERate",log=T,correlates.df)
{
  nse.fits <- list.dirs(model.fit.dir,recursive = F)
  names(nse.fits) <- c("1x","1000x","100x","10x")
  nse.dfs <- lapply(nse.fits,function(x)
  {
    params <- getParameterDataFrames(file.path(x,"final_restart"))
    nse.df <- params[[param.type]]
    hdis <- getHDI(file.path(x,"final_restart","R_objects","parameter.Rda"),param=param.type,log=log)
    nse.df <- nse.df %>% 
      left_join(hdis,by="Codon") %>%
      mutate(Quantile.Width = `97.5%` - `2.5%`,
             HDI.Width = HDI_high - HDI_low)
    return(nse.df)
  }) %>% 
    bind_rows(.id="Read.Depth")
  
  return(nse.dfs)
}


rescalePANSEWaitingTimes <- function(wait.time,
                                     target.mean.elong.rate = 9.3)
{
  ## Wait times are relative to methionine. rescale to codon/s
  elong.rate <- 1/wait.time$Wait.Time
  a <- (target.mean.elong.rate * sum(1/elong.rate))/61 # solve for a that will rescale elong rate such that harmonic mean matches target
  elong.rate.scaled <- a * elong.rate
  wait.time <- wait.time %>% 
    mutate(Wait.Time.Scaled = 1/elong.rate.scaled)
  
  return(wait.time)
}

calculateCDistribution <- function(wait.time, ## should contain wait times from PANSE rescaled to have units codon/sec as output by rescaleWaitingTimes
                                   deta.file = "../../Data/ROC_Estimates/saccharomyces_cerevisiae.max.cds_Selection.csv",
                                   Ne = 13600000,
                                   q = 4.19*10^-7)
{
  deta <- read_csv(deta.file) %>% 
    dplyr::rename(DEta = Mean)
  
  deta.w.wait.time <- deta %>% 
    left_join(wait.time,by=c("AA","Codon")) %>%
    mutate(Reference = ifelse(DEta == 0,1,0)) %>% 
    group_by(AA) %>% 
    mutate(Rel.WT = Wait.Time.Scaled - Wait.Time.Scaled[Reference == 1]) %>%
    filter(Reference != 1)
  
  dist.C <- deta.w.wait.time$DEta/(2*Ne*q*deta.w.wait.time$Rel.WT) 
  return(dist.C)
}

calculateMedianC <- function(dist.C,include.discrepancy.codons)
{
  if (include.discrepancy.codons)
  {
    C <- median(dist.C)
  } else {
    C <- median(dist.C[which(dist.C > 0)]) ## drop values that are negative, suggest disagreement between waiting times and \Delta\eta, use median
  }
  return(C)
}


#' A function for calculating the expected cost of translation including the direct cost of translation initiation and peptide elongation, ribosome pausing, and nonsense errors
#' @param panse.result.file.path A file path that points to a directory containing the results from a PANSE run
#' @param rfp A data.frame containing the rfp, similar to format used as input to PANSE
#' @param a.1 the initiation cost you want to use, Defaults to 4.
#' @param a.2 the direct elongation cost you want to use, Defaults to 4.
#' @param C The cost of pausing time in ATP/s C
#' @param w.0 The overhead indirect cost of translation initiation (units per sec)
#' @param target.mean.elong.rate The target harmonic mean of elongation rates

expectedCost <- function(panse.result.file.path,
                         rfp,
                         a.1 = 4,
                         a.2 = 4,
                         w.0= 10.8,
                         C = 5.5,
                         target.mean.elong.rate = 9.3)
{
  panse.result <- getParameterDataFrames(panse.result.file.path)
  
  wait.time <- panse.result$Wait.Time %>% 
    dplyr::rename(Wait.Time=Mean)
  
  wait.time <- rescalePANSEWaitingTimes(wait.time, 
                                        target.mean.elong.rate)
  rfp <- rfp %>%
    left_join(wait.time %>%
                dplyr::select(Codon,Wait.Time,Wait.Time.Scaled),
              by="Codon") %>%
    left_join(panse.result$NSERate %>% 
                dplyr::select(Codon,Mean) %>% 
                dplyr::rename(NSERate = Mean),
              by="Codon")
  rfp <- calculateSuccessProbability(rfp)
  
  rfp <- rfp %>%
    group_by(GeneID) %>%
    mutate(Sigma.prev.pos=dplyr::lag(Sigma.i,default=1),
           Wait.Time.Scaled.prev.pos = dplyr::lag(Wait.Time.Scaled,default=0),
           Total.Wait.time = cumsum(Wait.Time.Scaled.prev.pos)) %>%
    mutate(Pausing.cost.i = C * (w.0 + Total.Wait.time),
           Variable.cost.per.pos = (a.1 + (a.2 * (Position - 1)) + Pausing.cost.i) * Sigma.prev.pos * Pr.NSE,
           Variable.direct.cost.per.pos = (a.1 + (a.2 * (Position - 1))) * Sigma.prev.pos * Pr.NSE,
           Variable.indirect.cost.per.pos = (Pausing.cost.i) * Sigma.prev.pos * Pr.NSE
    )
  
  expected.cost <- rfp %>% 
    group_by(GeneID) %>%
    summarize(Length = n(),
              Sigma.n = prod(Prob.Successful.Per.Pos,na.rm=T),
              Expected.variable.cost.per.nse = sum(Variable.cost.per.pos,na.rm = T)/(1 - Sigma.n),
              Fixed.indirect.cost = C * (w.0 + sum(Wait.Time.Scaled)),
              Variable.direct.cost = (1/Sigma.n) * sum(Variable.direct.cost.per.pos),
              Variable.indirect.cost = (1/Sigma.n) * sum(Variable.indirect.cost.per.pos)) %>%
    mutate(Fixed.direct.cost = a.1 + a.2 * Length,
           Fixed.cost = Fixed.direct.cost + Fixed.indirect.cost,
           Variable.cost = Variable.direct.cost + Variable.indirect.cost,
           Direct.cost = Fixed.direct.cost + Variable.direct.cost,
           Indirect.cost = Fixed.indirect.cost + Variable.indirect.cost,
           Expected.cost.overall = Fixed.cost + Variable.cost
    ) %>%
    mutate(Result.file.path = panse.result.file.path,
           Initiation.cost.a1 = as.character(a.1), # convert to character for joining data.frames 
           Elongation.cost.a2 = as.character(a.2), 
           C = as.character(C),
           w.0 = as.character(w.0),
           Target.elongation.rate = as.character(target.mean.elong.rate))
  return(expected.cost)
}



getHDI <- function(parameter.path,param = "NSERate",log=T)
{
  parameter <- loadParameterObject(parameter.path)
  trace <- parameter$getTraceObject()
  if (param == "NSERate")
  {
    index <- 2
  } 
  if (param == "Alpha")
  {
    index <- 0
  }
  if (param == "Lambda")
  {
    index <- 1
  }
  param.trace <- trace$getCodonSpecificParameterTrace(index)
  codon.vec <- AnaCoDa::codons()[-(62:64)] # remove stop codons
  
  param.trace <- param.trace[[1]]
  names(param.trace) <- codon.vec
  param.hdi.95 <- lapply(param.trace,function(x){
    if(log)
    {
      hdi(log10(x))
    } else {
      hdi(x)
    }
  }) %>% 
    bind_rows(.id="Codon") %>%
    as_tibble() %>%
    dplyr::select(Codon,CI_low,CI_high) %>%
    dplyr::rename(HDI_low=CI_low,HDI_high=CI_high)
  median.value <- unlist(lapply(param.trace,function(x){
    if(log)
    {
      median(log10(x))
    } else {
      median(x)
    }
  }))
  
 param.hdi.95 <- param.hdi.95 %>%
    mutate(Median=median.value)
  return(param.hdi.95)
}

