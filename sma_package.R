###########################################################################
# Statistics for Microarray Analysis
# Exploratory analysis - Mainly preprocessing.
#
# Date : August 9, 2000
# Last update : May 17, 2001
#
# History:
#   May 17, 2001: Fix to norm.scale.func
#   March, 19: Splitting Rarray in to smaller files.  
#              Including Comments at the start of each function.
#   Nov, 20: Change the argument on plot.mva...it's not usable otherwise.
#            Bug fix ma.func
#   Nov, 13: Ben's Bug fix on stat.ma
#   Nov, 10: Change data structure from matrix to list of matrix.  
#   Sept, 28: Bug fix: ma.func
#   Feb 20, 2003 - bug fix to ma.func (As suggested by G. Smyth)
#   Apr 27, 2003 - fix bug in ma.func (when both R-Rb and G-Gb are negative should M give NA)
#
# Authors: Sandrine Dudoit and Yee Hwa (Jean) Yang.
# http://www.stat.berkeley.edu/~terry/zarray/Software/smacode.html
##########################################################################


##########################################################################
#  stat.gnames
#  History:  
#     March 19, 2001:  remove infinite values from the ordering.
#
##########################################################################

#########################################################################/**
# \name{stat.gnames}
# 
# \alias{stat.gnames}
# 
# \title{Sort Genes According to the Value of a Statistic}
# 
# \description{
# Lists genes and corresponding statistics in decreasing order of the
# statistics. This function applies to any type of statistic, including
# log ratios, one and two-sample t-statistics, and F-statistics. Missing
# values are ignored, as in \code{\link{sort}(..., na.last=NA)}. 
# }
# 
# \usage{
# stat.gnames(x, gnames, crit=0.05)
# }
# 
# \arguments{
#  \item{x}{a numeric vector containing the statistics for each
#  gene. Missing values (NAs) are allowed. }
#  
# \item{gnames}{a character vector containing the gene names.}
# 
#  \item{crit}{specifies the number of genes to be returned. If crit <
#  1, the crit*100\% genes with the largest x values are listed. If crit
#  >= 1, the crit genes with the largest x values are listed. }
# }
# 
# \value{
# List containing the following components 
#   \item{gnames}{gene names sorted in decreasing order of the
#  statistics in x.}
#  \item{t}{statistics sorted in decreasing order.}
# }
# 
# \author{
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} }
# 
# \seealso{\code{\link{stat.t2}}, \code{\link{order}}, \code{\link{sort}}.}
# 
# \examples{
# ## Calculating log ratio and performing a t test.
# data(MouseArray)
# ## mouse.setup <- init.grid()
# ## mouse.data <- init.data() ## see \emph{init.data}
# mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# cl <- c(rep(1,3), rep(2,3))
# mouse.t2 <- stat.t2(mouse.lratio, cl)
# 
# ## Looking at gene names
# ## Finding the top 10 t-statistics
# stat.gnames(abs(mouse.t2$t), mouse.gnames, crit=10)
# 
# ## Finding the top 1% of t-statistics
# stat.gnames(abs(mouse.t2$t), mouse.gnames, crit=0.01)
# 
# ## Finding the 10 extreme M values in the first slide
# stat.gnames(abs(mouse.lratio$M[, 1]), mouse.gnames, crit=10)
# }
# 
# \keyword{microarray.}
#*/#########################################################################

stat.gnames<-function(x, gnames, crit=0.05)
{
    ind <- is.infinite(x)
    x <- x[!ind]
    if (crit < 1) {
        which <- rev(order.na(x, na.last = FALSE))[1:(round(length(x) * 
            crit))]
        if (sum(is.na(x)) > (length(x) - round(length(x) * crit))) 
            warning("NA exists under your selection criteria")
    }
    if (crit >= 1) {
        which <- rev(order.na(x, na.last = FALSE))[1:crit]
        if (sum(is.na(x)) > (length(x) - crit)) 
            warning("NA exists under your selection criteria")
    }
    if (is.matrix(gnames) | is.data.frame(gnames)) 
      {
	gnames <- gnames[!ind, ]
        res <- list(gnames = gnames[which, ], t = x[which])
      }
    if (is.vector(gnames)) 
      {
	gnames <- gnames[!ind]
        res <- list(gnames = gnames[which], t = x[which])
      }
    res
}


##########################################################################
# Calculation of log ratios and normalization
# History:
#    March, 19, 2001: Modify the function so that stat.ma works with data 
#                     with no background values.
##########################################################################

#########################################################################/**
# \name{stat.ma}
# 
# \alias{stat.ma}
# \alias{ma.func}
# \alias{norm.l.func}
# \alias{norm.pin.func}
# \alias{norm.scale.func}
# 
# \title{Calculation of log Intensity Ratios and Average log Intensities}
# 
# \description{
# Computes the log intensity ratio \eqn{M = log_2 (R/G)} and the mean log
# intensity \eqn{A = log_2 \sqrt{RG}}{A = log_2(R*G)/2}, where R and G
# represent the fluorescence
# intensities in the red and green channels, respectively. Logarithms base
# 2 are used instead of natural or decimal logarithms as intensities are
# typically integers between 1 and \eqn{2^{16}}. The log intensity
# ratios M are normalized using one of the five available methods. 
# }
# 
# \usage{
# stat.ma(RG, layout, norm="p", pout=FALSE, ...)
# }
# 
# \arguments{
#   \item{RG}{
#     a list with 4 elements, each represents a matrix with p rows for p
#     genes and n columns for n slides. \cr
#     The first element "R" contains the raw red intensities from slide
#     i=1,...,n .\cr
#     Similarly, the second element "G" contains the raw green
#     intensities. \cr
#     The third element "Rb"  contains the background red intensities and \cr
#     the fourth element "Gb" contains the  background green intensities.\cr
#     This list structure can be generated by the interactive function
#     \code{\link{init.data}}. }
#   
#   \item{layout}{a list specifying the dimensions of the spot matrix
#   and the grid  
#     matrix.  This can be generated by calling \code{\link{init.grid}}.}
# 
#   \item{norm}{Character string, one of "n", "m", "l", "p" or "s".  This
#     argument specifies the type of normalization method to be performed:
#     "n" no normalization between the 2 channels; "m"
#     \code{\link{median}} normalization, which sets the median of log
#     intensity ratios to zero; "l" global \code{\link{lowess}}
#     normalization; "p" print-tip group lowess normalization and "s"
#     scaled print-tip group lowess normalization. The default method is
#     set to print-tip normalization.}
#
#   \item{pout}{if TRUE, an M vs. A plot will be produced. Otherwise,
#   a matrix of log intensity ratios and average log intensities is
#   return.  By default pout is set to FALSE.  The option pout='TURE'
#   is not yet implemented.}
#
#   \item{\dots}{other parameters used in \code{\link{ma.func}}. }
# }
# 
# \value{
#   List containing the following components:
#   
#   \item{M}{Matrix of log expression ratios \eqn{M = log_2 (R/G)}}
#   \item{A}{Matrix of average log intensities \eqn{A = log_2
#       \sqrt{RG}}{A = log_2(R*G)/2}}
#   For the matrix in each of the components, rows correspond to genes
#   and columns correspond to different hybridizations, that is
#   different slides.  
# }
# 
# \references{S. Dudoit, Y. H. Yang, M. J. Callow, and T. P. Speed. Statistical
# methods for identifying differentially expressed genes in replicated
# cDNA microarray experiments (Statistics, UC Berkeley, Tech Report \# 578).  }
# 
# \note{\code{\link{ma.func}}, \code{\link{norm.l.func}},
#   \code{\link{norm.scale.func}} and \code{\link{norm.pin.func}} are called by \code{\link{stat.ma}} and are not typically used on their own.}
# 
# \author{
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} \cr
#   Natalie Roberts, \email{nroberts@wehi.edu.au}
# }
# 
# \seealso{\code{\link{ma.func}}, \code{\link{norm.l.func}},
#   \code{\link{norm.pin.func}}, \code{\link{norm.scale.func}}, \code{\link{plot.mva}}, \code{\link{lowess}}.}
# 
# \examples{
# data(MouseArray)
# ## mouse.setup <- init.grid() 
# ## mouse.data <- init.data() ## see \emph{init.data} 
# mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# }
# 
# \keyword{microarray, log ratio.}
#*/#########################################################################

stat.ma <- function(RG, layout, norm="p", pout=FALSE, ...)
{
  n <- ncol(RG$R)
  res <- list(A=NULL, M=NULL)

  for(i in (1:n))
    {
##       RG <- apply(RG, 2, as.numeric)
      if(is.null(RG$R[,i])){
	stop(" Error: No data is given in RG$R \n")
      }
      if(is.null(RG$G[,i])){
	stop(" Error: No data is given in RG$G\n")
      }
      if(pout)
        stop("pout=TRUE is not implemented")
      tmp <-ma.func(R=RG$R[,i],G=RG$G[,i],Rb=RG$Rb[,i], Gb=RG$Gb[,i], layout=layout, norm=norm, pout=pout, ...)
      res$A<-cbind(res$A, tmp$A)
      res$M<-cbind(res$M, tmp$M)
    }
  res
}

##########################################################################
#  stat.norm.exp
#  History:  
#     March 19, 2001:  Original
#
##########################################################################

#########################################################################/**
# \name{stat.norm.exp}
# \alias{stat.norm.exp}
# \title{Normalization of log Intensity Ratios across slides / experiments.}
# 
# \description{
# Performs scale normalization across slides (experiments)}
# }
# 
# \usage{
# stat.norm.exp(X)
# }
# 
# \arguments{
#   \item{X}{X is a matrix of log intensity ratios \eqn{M=\log_2 (R/G)}
#   The rows of X correspond to genes and columns correspond to different 
#   hybridizations, that is different slides (experiments). 
# }
# 
# \value{
#   A matrix of normalized log intensity ratios across different slides. 
#   For the matrix in each of the components, rows correspond to genes
#   and columns correspond to different hybridizations, that is different 
#   slides.  This methods scale the matrix such that each column has the 
#   same median absolute deviation.
# }
# 
# \references{Y. H. Yang, S. Dudoit, P. Luu and T. P. Speed. 
#  Normalization for cDNA Microarray Data. (Statistics, UC Berkeley, 
#  Tech Report \# 589).  } 
# 
# \author{
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
# }
# 
# \seealso{\code{\link{ma.func}}, \code{\link{norm.l.func}},
#   \code{\link{norm.pin.func}}, \code{\link{norm.scale.func}}, \code{\link{plot.mva}}, \code{\link{lowess}}.}
# 
# \examples{
# data(MouseArray)
# ## mouse.setup <- init.grid() 
# ## mouse.data <- init.data() ## see \emph{init.data} 
# mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# mouse.norm.lratio <- stat.norm.exp(mouse.lratio$M)
# }
# 
# \keyword{microarray, log ratio, normalization.}
#*/#########################################################################

stat.norm.exp <- function(X)
  {
    n <- ncol(X)
    
    xmat.mad <- apply(X, 2, mad, na.rm=TRUE)
    
    denom <- (prod.na(xmat.mad))^(1/n)
    si <- xmat.mad / denom

    t(t(X) / si)
  }


##########################################################################
# Internal functions call by stat.ma
##########################################################################

##########################################################################
# ma.func
#
# March 20,  Remove pch="." from the code and set as an argument.
#            
##########################################################################                 
ma.func <- function (R, G, Rb, Gb, layout, norm = "p", pout, f = 0.3, extra.type="tci", crit1=0.025,crit2=crit1, nclass=10, labs=NULL, plot.type="n", col.ex=NULL, pch=pch, ...){
###
# extra.type ="t" for txt, ="p" for points, ="tci" for text ci, ="pci" for points ci ="lci" for lines confidence bands   
# crit is the size of pointwise confidence bands
# nclass (0 < nclass <1) is the proportion of points in each band i.e. smoothness of confbands
# plot.type ="n" plot normalised, ="r" raw data, ="b" both.
###                 
  if(is.null(Gb))
    cy3 <- G
  else
    cy3 <- G - Gb
  if(is.null(Rb))
    cy5 <- R
  else
    cy5 <- R - Rb

  A <- oA <- (log.na(cy3,2) +   log.na(cy5, 2))/2   #<- log.na(cy3 * cy5, 2)/2
  oM <- (log.na(cy5,2) - log.na(cy3,2))    #log.na(cy5/cy3, 2)
  if (norm == "n")
    M <- oM
  if (norm == "m")
    M <- oM - median(oM, na.rm = TRUE)
  if (norm == "l")
    M <- norm.l.func(oA, oM, f = f)
  if (norm == "p")
    M <- norm.pin.func(oA, oM, layout, f = f)
  if (norm =="s"){
    temp <- norm.pin.func(oA, oM, layout, f = f)
    M <- norm.scale.func(temp, layout)
  }
  if (pout){
    if(is.null(labs)) labs <- as.character(1:length(M))
    if(plot.type=="b") par(ask=TRUE) else par(ask=FALSE)
    if( ((plot.type == "b") | (plot.type == "r")) ){
##      par(mfrow = c(2, 1))
      plot(oA, oM, xlab = "A", ylab = "M", pch=pch, ...)
      if(extra.type== "t") text(oA,oM,labs,col=col.ex,...)
      if(extra.type== "p") points(oA,oM,col=col.ex,...)
      if(extra.type== "tci") plot.confband.text(oA,oM,crit1,crit2,nclass, labs,col=col.ex,...)
      if(extra.type== "pci") plot.confband.points(oA,oM,crit1,crit2,nclass,col=col.ex,...)
      if(extra.type== "lci") plot.confband.lines(oA,oM,crit1,crit2,nclass,col=col.ex,...)
      plot.smooth.line(oA, oM, f = 0.3)
    }
    if( ((plot.type == "b") | (plot.type =="n")) ){
      plot(A, M, xlab = "A", ylab = "Normalized M",pch=pch,...)
      if(extra.type== "t") txt(A,M,labs,col=col.ex,...)
      if(extra.type== "p") points(A,M,col=col.ex,...)
      if(extra.type== "tci") plot.confband.text(A,M,crit1,crit2,nclass,labs, col=col.ex,...)
      if(extra.type== "pci") plot.confband.points(A,M,crit1,crit2,nclass,col=col.ex,...)
      if(extra.type== "lci") plot.confband.lines(A,M,crit1,crit2,nclass,col=col.ex,...)
      plot.smooth.line(A, M, f = 0.3)
    }
    par(ask=FALSE)
##    par(mfrow = c(1, 1))
  }
  else{ list(M = M, A = A) }
}

norm.scale.func <- function(x, layout, x.names=NULL)
  {
    n <- layout$nspot.r * layout$nspot.c * layout$ngrid.r * layout$ngrid.c
    nperpin <- layout$nspot.r * layout$nspot.c
    npin <- layout$ngrid.r * layout$ngrid.c
    
    if(is.vector(x)){
      if((length(x) != n) & (is.null(x.names)))
        {
stop(" Error: Length of vector different from total number of spots and vector has no row.name.\n")
        }
      if ((length(as.vector(x)) != n) & (!is.null(x.names)))
        {
          y <- x; x <- rep(NA, n);
          x[as.integer(x.names)] <- y
        }
      xmat <- matrix(x, nrow = nperpin)
      vect <- TRUE
         } 

    if(is.matrix(x)){
      xmat <- x
      vect <- FALSE
                    }
    
    xmat.mad <- apply(xmat, 2, mad, na.rm=TRUE)
    sigma2 <- (1/nperpin) * exp((1/npin)*sum.na(log(xmat.mad)))
    si <- xmat.mad / (sigma2 * nperpin)
    xmat.s <- t(t(xmat) / si)

  if(vect)
      res <- as.vector(xmat.s)
    else
      res <- xmat.s
    res
  }

norm.l.func <- function(A, M, ...)
  {
    ind <- is.na(A) | is.na(M) | is.infinite(A) | is.infinite(M)
    smoothnum <- lowess(A[!ind], M[!ind], ...)
    lowesslratio <- M
    lowesslratio[!ind] <- M[!ind] - approx(smoothnum, xout = A[!ind])$y
    lowesslratio
  }

norm.pin.func <- function(A, M, layout, ...)
  {
    npin <- layout$ngrid.r * layout$ngrid.c
    pin <- c(0, rep(layout$nspot.r * layout$nspot.c, npin) * (1:npin))
    ind <- 1:length(M)
    lowessratio <- M
    for(j in 1:npin) {
      index <- ((pin[j] + 1) <= ind) & (ind <= pin[j + 1])
      tM <- M[index]
      tA <- A[index]
      ind2 <- is.na(tM) | is.na(tA) | is.infinite(tM) | is.infinite(tA)
      smoothnum <- lowess(tA[!ind2], tM[!ind2], ...)
      lowessratio[index][!ind2] <- tM[!ind2] - approx(smoothnum, xout = tA[!ind2])$y
    }
    lowessratio
  }


##########################################################################
#                                End of file
##########################################################################
###########################################################################
# Statistics for Microarray Analysis
# Bayesian Method
#
# Date : October 2, 2001
#
# History:
#
# Authors: Ingrid Lönnstedt  and Yee Hwa (Jean) Yang.
# Written by Ingrid with help from Jean.
##########################################################################

#stat.bayesian calculates a lodscore (lods) for each gene in an experiment, using the normalized M-values (output from stat.ma), the number of slides (nb), and the number of replicates for each gene within each slide (nw). If there are j replicates within slides, the vectors of M-values for each slide should be on the form M11, ..., M1j, M21, ...M2j, ..., Mgj, where g is the number of genes.

#stat.bay.est calculates a lodscore (lods) for each gene in an experiment, using independent, sufficient statistics for the effect and its variance of each gene. 

#plot.bayesian plots the results of any of the above, highlighting genes which meet userdefined criteria.

######################################################
# Main Program
######################################################

stat.bayesian <- function(M=NULL, nb=NULL, nw=1, Xprep=NULL, para=list(p=0.01, v = NULL, a=NULL, c = NULL,k=NULL))
  {
    ## Input M (output from stat.ma) as well as nb (and nw if nw>1)
    ## Xprep and para are calculated and used for calculating lods.
    ## If Xprep is given in the function input, M, nb and nw are unnecessary.

    ## Setting up Data
    if(is.null(Xprep))
      {
    	Xprep<- setup.bayesian(M=M, nb=nb, nw=nw)
      }
    nb<-Xprep$nb
    nw<-Xprep$nw
    Mbar<-Xprep$Mbar
    SSW<-Xprep$SSW
    SSB<-Xprep$SSB

    ## Setting up parameters
     if(is.null(para$v) | is.null(para$a)) 
     {
	va <- va.func(Vest=SSB/(nb-1), k=nb)
	if (is.null(para$v)) para$v<-va$v
	if (is.null(para$a)) para$a<-va$a
     }

     if(is.null(para$k)) 
     {
	if(is.null(SSW)) para$k<-0  else para$k<-median(SSB/(nb-1)/SSW*nb*(nw-1))
     }

     if(is.null(para$c)) para$c<-c.func(Xprep=Xprep, para=para)

     lods<-lods.func(Xprep, para)
     list(Xprep=Xprep, lods=lods,para=para)
}


stat.bay.est <- function(M=NULL, Xprep=NULL, para=list(p=0.01, v = NULL, a=NULL, c = NULL))
  {
    ## Input Xprep or M (output from stat.ma).
    ## If M is given, stat.bay.est assumes the experiment consists of ncol(M) microarray slides all measuring the same effect (which will be stimated by Mbar)
    ## Para is calculated and used for calculating lods.
    
    ## Setting up Data
    if(is.null(Xprep))
      {
    	Xprep$Mbar<-apply(M,1,mean.na)   #Effect estimate
        Xprep$Vest<-apply(M,1,var.na)    #Variance estimate
        Xprep$k<-ncol(M)                 #Variance constant
        Xprep$f<-Xprep$k-1               #Degrees of freedom
      }
    Mbar<-Xprep$Mbar
    Vest<-Xprep$Vest
    k<-Xprep$k
    f<-Xprep$f

    ## Setting up parameters
     if(is.null(para$v) | is.null(para$a)) 
     {
	va <- va.func(Vest=Vest, k=k)
	if (is.null(para$v)) para$v<-va$v
	if (is.null(para$a)) para$a<-va$a
     }

     if(is.null(para$c)) para$c<-c.func.est(Xprep=Xprep, para=para)

     lods<-lods.func.est(Xprep, para)
     list(Xprep=Xprep, lods=lods,para=para)
}


plot.bayesian<-function(x=NULL,Mbar=x$Xprep$Mbar,lods=x$lods, type="t",spec=50, ch=NULL, col='black',...){
  #bay <- x
  index<-NULL
  if(type=="t") index<-(1:length(lods))[lods>sort(lods[!is.na(lods)])[length(lods[!is.na(lods)])-spec]]   
  if(type=="c") index<-(1:length(lods))[lods >= spec]
  if(type=="i") index<-spec                       
  plot(Mbar, lods, xlab="Effect estimate", ylab="Lodsratio", main="Lodsratio vs Effect estimate", type="n")
  if (!is.null(index)){
  points(Mbar[-index], lods[-index], pch='.')    
  if(is.null(ch))  text(Mbar[index],lods[index],labels=index,col=col)  else points(Mbar[index],lods[index],pch=ch,col=col)
}
  if(is.null(index)){points(Mbar, lods, pch='.')}
}
  
    
######################################################
#Functions
######################################################

setup.bayesian <- function(M, nb=NULL, nw=1)
  {
    if (nw == 1)
    {
	nb<-ncol(M)
	SSW<-NULL
	SSB<-apply(M,1,var.na)*(nb-1)
	Mbar<-apply(M,1,mean.na)
    }

    if (nw > 1)
    {
	
	if (nb > 1)
	{
  	  Mtmp<-NULL
	  for (i in 1:nb)
	  {
	    for (j in 1:nw)
	    {
	      Mtmp<-cbind(Mtmp,M[seq(j,nrow(M),nw),i])
	    }
	  }
	  Mbar<-apply(Mtmp,1,mean.na)

	  Mslide<-NULL
	  SSW<-rep(0,nrow(M)/nw)
	  SSB<-rep(0,nrow(M)/nw)
	  for (i in 1:nb)
	  {
	    Mslide<-cbind(Mslide,apply(Mtmp[,((i-1)*nw+1):(i*nw)],1,mean.na))
	    SSW<-SSW+apply((Mtmp[,((i-1)*nw+1):(i*nw)]-Mslide[,i])**2,1,sum.na)	
	    SSB<-SSB+nw*((Mslide[,i]-Mbar)**2)
	  }
	}

	if (nb == 1)
	{

  	  Mtmp<-NULL

          for (j in 1:nw)
          {
	    Mtmp<-cbind(Mtmp,M[seq(j,length(M),nw)])
	  }

	  SSB<-apply(Mtmp,1,var.na)*(nw-1)
	  nb<-nw
	  nw<-1
	  SSW<-NULL
	  Mbar<-apply(Mtmp,1,mean.na)
	}
	
    }	
	list(Mbar=Mbar, SSB=SSB, SSW=SSW, nb=nb, nw=nw)
  }

####################################################################
#Lods (or lods.func.est) calculates the logodds ratio.

lods.func<-function(Xprep=list(Mbar=Mbar, SSB=SSB, SSW=SSW, nb=nb, nw=nw), para=list(p=p, c=c, v=v, a=a, k=k)){
   Mbar <- Xprep$Mbar            #overall means
   nb<-Xprep$nb
   nw<-Xprep$nw
   SSB <- Xprep$SSB              #sums of squares between slides
   SSW <- Xprep$SSW              #sums of squares within slides  
   if(is.null(Xprep$SSW)) SSW <- rep(0, length(SSB))

   p <- para$p
   v <- para$v
   a <- para$a
   c <- para$c
   k <- para$k

   odds1<-p/(1-p)
   odds2<-1/(1+nb*nw*c)
   odds3<-a + (SSB + k*SSW)/nw/nb   
   odds4<-(odds3+Mbar^2)/(odds3+Mbar^2/(1+nb*nw*c))
   odds<-odds1*(odds2**(1/2))*(odds4**(v+nw*nb/2))
   log(odds)
 }

lods.func.est<-function(Xprep=list(Mbar=Mbar, Vest=Vest, k=k, f=f), para=list(p=p, c=c, v=v, a=a)){

   Mbar<-Xprep$Mbar
   Vest<-Xprep$Vest
   k<-Xprep$k
   f<-Xprep$f

   p <- para$p
   v <- para$v
   a <- para$a
   c <- para$c

   odds1<-p/(1-p)
   odds2<-1/(1+k*c)
   odds3<-a + f*Vest/k  
   odds4<-(odds3+Mbar^2)/(odds3+Mbar^2/(1+k*c))
   odds<-odds1*(odds2**(1/2))*(odds4**(v+f/2+1/2))
   log(odds)
 }

######################################################

#va.func estimates v and a by the method of moments, so that a*k/(2*sigma^2) ~Gamma(v,1), Vest are the genewise estimates of sigma^2 and k a constant such that the expected variance of the effect estimates are sigma^2/k. 

va.func<-function(Vest, k){
av.var<-mean.na(Vest)
var.var<-sum.na((Vest-mean.na(Vest))^2)/(length.na(Vest)-1)
vhat<-(2*var.var+av.var^2)/var.var
ahat<-av.var/k*2*(vhat-1)
list (v=vhat, a=ahat)
}


######################################################
#c.func (or c.func.est) uses ls.variance and sq.func to estimate c. It uses a least squares estimate so that for all Mbar-values, Mbar~N(0,simga^2/k) and for the top p proportion of the genes, the averages are ~N(0,c*sigma^2).

c.func<-function(Xprep, para){

#Estimate the variance sigma²
  var.est<-mean.na(Xprep$SSB/(Xprep$nb-1))
  start<-var.est/10/Xprep$nb
  end<-var.est*10/Xprep$nb
  
  sigma2<-ls.variance(X=Xprep$Mbar[!is.na(Xprep$Mbar)], var.start=start, var.stop=end, nclass=100)*Xprep$nb
  
  sigma2


#Estimate c*sigma²
  para$c<-1.5
  if (is.null(Xprep$SSW)) Xprep$SSW<-rep(0,length(Xprep$SSB))
  l<-stat.bayesian(Xprep=Xprep, para=para)$lods
  top.set<-(1:(length(Xprep$Mbar)/Xprep$nw))[!is.na(l)][rank(l[!is.na(l)])>(length.na(l)-round(length.na(l)*para$p))]

  var.est<-var.na(Xprep$Mbar[top.set])
  start<-var.est/10
  end<-var.est*10
  
  csigma2<-ls.variance(X=Xprep$Mbar[top.set], var.start=start, var.stop=end, zeros=FALSE)
  
  csigma2/sigma2
}


c.func.est<-function(Xprep, para){

#Estimate the variance sigma²
  var.est<-mean.na(Xprep$Vest)
  start<-var.est/10/Xprep$k
  end<-var.est*10/Xprep$k
  
  sigma2<-ls.variance(X=Xprep$Mbar[!is.na(Xprep$Mbar)], var.start=start, var.stop=end, nclass=100)*Xprep$k
  
  sigma2


#Estimate c*sigma²
  para$c<-1.5
  l<-stat.bay.est(Xprep=Xprep, para=para)$lods
  top.set<-(1:(length(Xprep$Mbar)))[!is.na(l)][rank(l[!is.na(l)])>(length.na(l)-round(length.na(l)*para$p))]

  var.est<-mean.na(Xprep$Vest[top.set])
  start<-var.est/10
  end<-var.est*10
  
  csigma2<-ls.variance(X=Xprep$Mbar[top.set], var.start=start, var.stop=end, zeros=FALSE)
  
  csigma2/sigma2
}


ls.variance<-function(X, var.start, var.stop, var.steps=100, nclass=NULL, zeros=TRUE) {

  var.seq<-seq(var.start, var.stop, (var.stop-var.start)/var.steps)
  sq<-rep(0,length(var.seq))
  i<-0
  for (v in var.seq){
    i<-i+1
      sq[i] <- sq.func(X=X, var=v, nclass=nclass, zeros=zeros)
  }  
  
  minsq<-min(sq)
  for (i in 1:length(var.seq)){
      if (sq[i]==minsq)  min.v<-i
  }
  var.seq[min.v]
}

sq.func<-function(X,var,nclass=NULL, zeros=TRUE){

    if (is.null(nclass)) index<-1:length(X) else index <- seq(1, length(X), round(length(X)/nclass))
    hst<-hist(X[index], plot=FALSE, freq=FALSE, nclass=20)

    if (zeros) sumindex<-(1:length(hst$counts)) else sumindex<-(1:length(hst$counts))[hst$density > 0]

    theo<-dnorm(hst$mids[sumindex], mean=0 ,sd=sqrt(var))
    sum.na((hst$density[sumindex]-theo)**2)
}












###########################################################################
# Statistics for Microarray Analysis for R
# Discriminant analysis
#
# Date : August 21, 2000
# Last update : April 13, 2001
#
# Authors: Sandrine Dudoit, Yee Hwa (Jean) Yang, and Jane Fridlyand.
##########################################################################

##########################################################################
#                       A Red-Green Color Map
##########################################################################

########################################################################/**
#                            
# \name{rgcolors.func}
# 
# \alias{rgcolors.func}
# 
# \title{Red and Green Color Specification}
# 
# \description{
# This function creates a vector of n ``contiguous'' colors,
# corresponding to n intensities (between 0 and 1) of the red, green
# and blue primaries, with the blue intensities set to zero. The
# values returned by \code{rgcolors.func} can be used with a
# \code{col=} specification in graphics functions or in
# \code{\link{par}}.  
# }
# 
# \usage{
# rgcolors.func(n=50)
# }
# 
# \arguments{
#  \item{n}{the number of colors (>= 1) to be used in the red and
#  green palette. } 
# 
# }
# \value{a character vector of color names. Colors are specified
# directly in terms of their RGB components with a string of the form
# "\#RRGGBB", where each of the pairs RR, GG, BB consist of two
# hexadecimal digits giving a value in the range 00 to FF. 
#  }
# 
# 
# \author{
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} \cr
#   Jane Fridlyand, \email{janef@stat.berkeley.edu}
# }
# 
# \seealso{\code{\link{plot.cor}}, \code{\link{plot.mat}},
# \code{\link{colors}}, \code{\link{rgb}}, \code{\link{image}}.} 
# 
# \examples{
# rgcolors.func(n=5)
# ## The following vector is returned:
# ## "#00FF00" "#40BF00" "#808000" "#BF4000" "#FF0000"
# }
# 
# \keyword{Microarray, RGB image.}
# 
#*/#######################################################################
                          
rgcolors.func<-function(n = 50) 
{
  k <- round(n/2)     
  r <- c(rep(0, k), seq(0, 1, length = k))     
  g <- c(rev(seq(0, 1, length = k)), rep(0, k))     
  res <- rgb(r, g, rep(0, 2 * k))     
  res 
}               

##########################################################################
#                Images of data matrices and correlation matrices
##########################################################################
########################################################################/**
# \name{plot.cor}
# 
# \alias{plot.cor}
# 
# \title{Red and Green Color Image of Correlation Matrix}
# 
# \description{
# This function produces a red and green color image of a correlation
# matrix using an RGB color specification. Increasingly positive
# correlations are represented with reds of increasing intensity, and
# increasingly negative correlations are represented with greens of
# increasing intensity.  
# }
# 
# \usage{
# plot.cor(X, new=F, nrgcols=50, labels=FALSE, labcols=1, title="", ...)
# }
# 
# \arguments{
#  \item{X}{a matrix of numerical values.}
#  \item{new}{If \code{new=F}, \code{X} must already be a correlation
#  matrix. If \code{new=T}, the correlation matrix for the columns of
#  \code{X} is computed and displayed in the image.} 
#  \item{nrgcols}{the number of colors (>= 1) to be used in the red
#  and green palette.} 
#  \item{labels}{vector of character strings to be placed at the
#  tickpoints, labels for the columns of \code{X}.} 
#  \item{labcols}{colors to be used for the labels of the columns of
#  \code{X}. \code{labcols} can have either length 1, in which case
#  all the labels are displayed using the same color, or the same
#  length as \code{labels}, in which case a color is specified for the
#  label of each column of \code{X}.} 
#  \item{title}{character string, overall title for the plot.}
#  \item{\dots}{graphical parameters may also be supplied as arguments to
#           the function (see \code{\link{par}}). For comparison purposes, 
#  it is good to set \code{zlim=c(-1,1)}.}
# }
# }
# 
# 
# \author{
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu}
# }
# 
# \seealso{\code{\link{plot.mat}},\code{\link{rgcolors.func}},
# \code{\link{cor.na}}, \code{\link{cor}}, \code{\link{image}},
# \code{\link{rgb}}.} 
# 
# 
# \keyword{Microarray, correlation matrix, image.}
# 
# 
#*/#######################################################################

 plot.cor<-function(x, new=FALSE, nrgcols=50, labels=FALSE, labcols=1, title="", ...)
 {
#   X <- x
   n<-ncol(x)
   corr<-x
 
   if(new)
     corr<-cor.na(x)
  
   image(1:n,1:n,corr[,n:1],col=rgcolors.func(nrgcols),axes=FALSE, xlab="", ylab="",... ) 
 
  if(length(labcols)==1){
    axis(2,at=n:1,labels=labels,las=2,cex.axis=0.6,col.axis=labcols)
    axis(3,at=1:n,labels=labels,las=2,cex.axis=0.6,col.axis=labcols)
  }

  if(length(labcols)==n){
    cols<-unique(labcols)
    for(i in 1:length(cols)){
      which<-(1:n)[labcols==cols[i]]
      axis(2,at=(n:1)[which],labels=labels[which],las=2,cex.axis=0.6,col.axis=cols[i])
      axis(3,at=which,labels=labels[which],las=2,cex.axis=0.6,col.axis=cols[i])
     }
  }

  mtext(title,side=3,line=3)
  box()
}

########################################################################/**
# \name{plot.mat}
# 
# \alias{plot.mat}
# 
# \title{Red and Green Color Image of Data Matrix}
# 
# \description{This function produces a red and green color image of a
# data matrix using an RGB color specification. Larger entries are
# represented with reds of increasing intensity, and smaller entries
# are represented with greens of increasing intensity.  
# }
# 
# \usage{
# plot.mat(X, nrgcols=50, rlabels=FALSE, clabels=FALSE, rcols=1, ccols=1, title="",...)
# }
# 
# %- maybe also `usage' for other objects documented here.
# 
# \arguments{
#  \item{X}{a matrix of numbers.}
#  \item{nrgcols}{the number of colors (>= 1) to be used in the red
#  and green palette.} 
#  \item{rlabels}{vector of character strings to be placed at the row
#  tickpoints, labels for the rows of \code{X}.} 
#  \item{clabels}{vector of character strings to be placed at the
#  column tickpoints, labels for the columns of \code{X}.} 
#  \item{rcols}{colors to be used for the labels of the rows of
#  \code{X}. \code{rcols} can have either length 1, in which case
#  all the labels are displayed using the same color, or the same
#  length as \code{rlabels}, in which case a color is specified for the
#  label of each row of \code{X}.} 
#  \item{ccols}{colors to be used for the labels of the columns of
#  \code{X}. \code{ccols} can have either length 1, in which case
#  all the labels are displayed using the same color, or the same
#  length as \code{clabels}, in which case a color is specified for the
#  label of each column of \code{X}.} 
#  \item{title}{character string, overall title for the plot.}
#  \item{\dots}{graphical parameters may also be supplied as arguments  to
#           the function (see \code{\link{par}}).  E.g. \code{zlim=c(-3,3)}}
# }
# 
# %\references{ ~put references to the literature/web site here ~ }
# 
# 
# \author{
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu}
# }
# 
# \seealso{\code{\link{plot.cor}}, \code{\link{rgcolors.func}},
# \code{\link{cor.na}}, \code{\link{cor}}, \code{\link{image}},
# \code{\link{rgb}}.} 
# 
# \examples{
# data(MouseArray)
# ##mouse.setup <- init.grid()
# ##mouse.data <- init.data() ## see \emph{init.data}
# mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# 
# ## Looking at log ratios of mouse1
# plot.mat(spatial.func(mouse.lratio$M[,1], mouse.setup))
# }
# 
# \keyword{Microarray, image of data matrix.} 
# 
# 
#*/#######################################################################

plot.mat<-function(x, nrgcols=50, rlabels=FALSE, clabels=FALSE, rcols=1, ccols=1, title="", ...)
{
#  X <-x
  n<-nrow(x)
  p<-ncol(x)	  

  image(1:p,1:n,t(x[n:1,]),col=rgcolors.func(nrgcols),axes=FALSE, xlab="", ylab="", ... ) 

  if(length(ccols)==1){
    axis(3,at=1:p,labels=clabels,las=2,cex.axis=0.6,col.axis=ccols)
      }

  if(length(ccols)==p){
    cols<-unique(ccols)
    for(i in 1:length(cols)){
      which<-(1:p)[ccols==cols[i]]
      axis(3,at=which,labels=clabels[which],las=2,cex.axis=0.6,col.axis=cols[i])
     }
  }

  if(length(rcols)==1){
    axis(2,at=n:1,labels=rlabels,las=2,cex.axis=0.6,col.axis=rcols)
      }

  if(length(rcols)==n){
    cols<-unique(rcols)
    for(i in 1:length(cols)){
      which<-(1:n)[rcols==cols[i]]
      axis(2,at=(n:1)[which],labels=rlabels[which],las=2,cex.axis=0.6,col.axis=cols[i])
     }
  }

  mtext(title,side=3,line=3)
  box()
}


###########################################################################
###########################################################################
##                Functions for discriminant analysis
###########################################################################
##########################################################################

########################################################################/**
# \name{stat.bwss}
# 
# \alias{stat.bwss}
# 
# 
# \title{Between and Within Group Sum of Squares Calculation}
# 
# \description{This function computes the between and within group sum
# of squares for each row of a matrix which may have NAs. 
# }
# 
# \usage{
# stat.bwss(x, cl)
# }
# 
# \arguments{
#  \item{x}{a matrix, NAs allowed. In the microarray context, the rows
#  of \code{X} could correspond to genes and the columns to different
#  hybridizations.} 
#  \item{cl}{column labels, must be consecutive integers.} 
# }
# 
# \value{
# List containing the following components
#   \item{wn}{matrix with class sizes for each row of \code{X};}
#   \item{BW}{vector of BSS/WSS for each row of \code{X};}
#   \item{BSS}{between group sum of squares for each row of \code{X};}
#   \item{WSS}{within group sum of squares for each row of \code{X};}
#   \item{TSS}{total sum of squares for each row of \code{X};}
#   \item{tvar}{variance for each row of \code{X}.
# }
# 
# \references{S. Dudoit, J. Fridlyand, and T. P. Speed. Comparison of
# Discrimination Methods for the Classification of Tumors Using Gene
# Expression Data. June 2000. (Statistics, UC Berkeley,  Tech Report \#
# 576). } 
# 
# \author{
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu}  \cr
#   Jane Fridlyand, \email{janef@stat.berkeley.edu}
# }
# 
# \seealso{\code{\link{var.na}}, \code{\link{var}}, \code{\link{apply}}.}
# 
# \keyword{Microarray, BSS, WSS, sum of squares.}
#
#*/#########################################################################

 
 stat.bwss<-function(x, cl)
 {
 # Compute BSS/WSS for each row of a matrix which may have NA
 # Columns have labels cl=consecutive integers
 
   K <- max(cl) - min(cl) + 1
   tvar <- apply(x, 1, var.na)
   tn <- apply(!is.na(x), 1, sum)
   wvar <- matrix(0, nrow(x), K)
   wn <- matrix(0, nrow(x), K)
 
   for(i in (1:K)) {
     if(sum(cl == (i + min(cl) - 1)) == 1) {
       wvar[, i] <- 0
       wn[, i] <- 1
     }
     if(sum(cl == (i + min(cl) - 1)) > 1) {
       wvar[, i] <- apply(x[, cl == (i + min(cl) - 1)], 1, var.na)
       wn[, i] <- apply(!is.na(x[, cl == (i + min(cl) - 1)]), 1, sum)
     }
   }
 
   WSS <- apply(wvar * (wn - 1), 1, sum.na)
   TSS <- tvar * (tn - 1)
   BSS<- TSS - WSS
   BW <- BSS/WSS
 
   list(wn=wn,bw=BW,bss=BSS,wss=WSS,tss=TSS,tvar=tvar)
 }
 
########################################################################/**
# \name{stat.diag.da}
# 
# \alias{stat.diag.da}
# 
# \title{Diagonal Discriminant Analysis}
# 
# \description{
# This function implements a simple Gaussian maximum likelihood
# discriminant rule, for diagonal class covariance matrices.  
# }
# 
# \usage{
# stat.diag.da(ls, cll, ts, pool=1)
# }
# 
# \arguments{
#  \item{ls}{learning set data matrix, with rows corresponding to
# cases (i.e., mRNA samples) and columns to predictor variables 
# (i.e., genes).} 
#  \item{cll}{class labels for learning set, must be consecutive integers.}
#  \item{ts}{test set data matrix, with rows corresponding to cases
# and columns to predictor variables.} 
#  \item{pool}{logical flag. If \code{pool=1}, the covariance matrices
# are assumed to be constant across classes and the discriminant rule
# is linear in the data. If \code{pool=0}, the covariance matrices may
# vary across classes and the discriminant rule is quadratic in the
# data.} 
# } 
# 
# \value{
#  List containing the following components 
#   \item{pred}{vector of class predictions for the test set.}
# }
# 
# \references{
# S. Dudoit, J. Fridlyand, and T. P. Speed. Comparison of
# Discrimination Methods for the Classification of Tumors Using Gene 
# Expression Data. June 2000. (Statistics, UC Berkeley, Tech Report \#576).}
# 
# \author{
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} \cr
#   Jane Fridlyand, \email{janef@stat.berkeley.edu}
# }
# 
# %\seealso{ ~~objects to SEE ALSO as \code{\link{~~fun~~}}, ~~~ }
# 
# \keyword{Microarray, discriminant analysis, maximum likelihood, classification.}
# 
# 
#*/#########################################################################

stat.diag.da<-function(ls, cll, ts, pool=1)
{
  ls <- as.matrix(ls)
  ts <- as.matrix(ts)
  n <- nrow(ls)
  p <- ncol(ls)	

  nk<-rep(0,max(cll)-min(cll)+1)
  K<-length(nk)
  m<-matrix(0,K,p)
  v<-matrix(0,K,p)
  disc<-matrix(0,nrow(ts),K)

  # Class means and variances
  for(k in (1:K)) 
    {
      which <- (cll == k + min(cll) - 1)
      nk[k]<-sum.na(which)
      m[k,]<-apply(ls[which,  ], 2, mean.na)
      v[k,]<-apply(ls[which,  ], 2, var.na)
    }

  # Pooled estimates of variances
  vp<-apply(v,2,function(z) sum.na((nk-1)*z)/(n-K))
 
  if(pool==1) # LDA 
    {
      for(k in (1:K))
        disc[,k]<-apply(ts,1,function(z) sum.na((z-m[k,])^2/vp))
    }
  if(pool==0) # QDA
    {
      for(k in (1:K))
        disc[,k]<-apply(ts,1,function(z) (sum.na((z-m[k,])^2/v[k,]) + sum.na(log(v[k,]))) )
    }

  # predictions
  pred<-apply(disc,1,function(z) (min(cll):max(cll))[order.na(z)[1]])
  
  list(pred=pred)
}
 
###########################################################################
# Statistical Microarray Analysis for R 
# Exploratory analysis (i)
#      Initialization functions
#
# Date : August 9, 2000
# Last update : November 12, 2000
#
# History:
#   March, 19: Insert comments from the help files
#   Nov, 10: Change data structure from matrix to list of matrix.  
#   Feb 12, 2003: Fix a bug in init.name.exp  init.readexp -> init.read.exp
#   Aug 15, 2003: Allow teh column names to be specified in read.genepix
#                 because some newer versions of genepix have more
#                 columns. We will allow the user to decide how they want to
#                 name the columns.
#
# Authors: Sandrine Dudoit, Yee Hwa (Jean) Yang and Natalie Roberts
#          with occasional maintanence from B. M. Bolstad
##########################################################################

##########################################################################
# Read in data from Spot output file
##########################################################################

########################################################################/**
# \name{write.spot}
# 
# \alias{write.spot}
# 
# \title{Writing in Data Generated by the Image Analysis Package Spot}
# \description{
#   Function writes in a data file in a tab delimited table format.
#   }
# 
# \usage{
# write.spot(x, imageid, batch="output")
# }
# 
# \arguments{
#  \item{x}{the object to be written, typically a data frame. If not, it
#    is attempted to create one from it.}
#  \item{imageid}{integer value; the index of the slide which is
#    considered}
#  \item{batch}{character string, this refers to the name of a collection
#    of experiments. The default batch name is "output".}
# }
# \details{
#   This function writes the data in for each imageid, assigning each file
#   the filename which takes the default form of "output".imageid.spot. The
#   column names of x are written along with x in the table format
# }
# 
# \references{ Spot manual \url{
#   http://www.cmis.csiro.au/iap/Spot/spotmanual.htm}}  
#   }
# 
# \author{ Jessica Mar }
# 
# \seealso{ \code{\link{write.table}, \code{\link{init.read.spot}} }
# 
# \examples{## Setting up the data
# ## library(Spot)
# ## SetParameters("mouse")
# ## Here is what you should see:
# ## Enter number of rows of grids per image (ngrid.r): 4
# ## Enter number of columns of grids per image (ngrid.c): 4
# ## Enter number of rows of spots per grid (nspot.r): 19
# ## Enter number of columns of spots per grid (nspot.c): 21
# ## Enter top/bottom translation tolerance, default is 50 (tolerance.r): 20
# ## Enter left/right translation, default is 50 (tolerance.c): 30
# ## Initialization complete
# 
# ##Inputting Image Data 
# ## SetImages("mouse")
# ## Combining the red and green channels for the first slide
# ## mouse.array <- Spots("mouse", 1)
# 
# ## Calling the function to write the data in
# ## write.spot(mouse.array, 1, "mouse")
# }
# 
# \keyword{microarray, Spot, Genepix.} 
# 
#*/#####################################################################
 
write.spot <- function(x, imageid, batch="output")
  {
    if(is.numeric(imageid))
      {
        newname <- paste(batch, imageid, "spot", sep=".")
      }
    if(is.character(imageid))
      newname <- imageid
    if(!is.character(imageid)&!is.numeric(imageid))
      {
        stop("Warning: imageid must be a number or a character")
      } 
    write.table(x, newname, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
  }

########################################################################/**
# \name{read.spot}
# 
# \alias{read.spot}
# 
# \title{Reading in Data Generated by the Image Analysis Package Spot}
# 
# \description{
# Reads in a data file in table format and creates a data frame with the
# same number of rows as there are lines in the file, and the same number
# of columns as there are fields in the file.\cr
# `read.spot': reads in the data file generated by the microarray image
#  extraction library Spot. 
# }
# 
# \usage{
# read.spot(name, dir=".", sep=",", header=T, ...)
# }
# 
# \arguments{
#  \item{name}{character string naming the data file from which to read the
# data. }
# 
#  \item{dir}{character string naming the directory that contains the data
#    file. The default setting is the current directory.}
# 
#  \item{sep}{the field separator (single character), often "\t"  for
#  tab delimited fields.  If omitted,  any amount of white space
#  (blanks or tabs) can separate the fields.  To read fixed format
#  files, make  sep  a numeric vector giving the initial columns of
#  the fields. } 
# 
#  \item{header}{logical flag: if TRUE, then the first line of the
#  file is used as the variable names of the resulting data frame. }
# 
#  \item{\dots}{parameters for read.table may also be supplied as arguments to
#  the function (see  \code{\link{read.table}}).  }
# }
# 
# \value{as in \code{\link{read.table}}, a data frame
# (\code{\link{data.frame}}) containing a representation of the data
# in the file.  
# }
# 
# \seealso{\code{\link{read.table}}, \code{\link{data.frame}},
#   \code{\link{write.spot}}, \code{\link{read.genepix}}.}
# 
# \examples{

# ## write.spot(mouse.array, 1, "mouse")
# ## mouse1 <- read.spot("mouse.1.spot")
# }
# 
# \keyword{microarray, Spot, GenePix.}
# 
#*/#####################################################################

read.spot <- function(name, dir=".", sep="\t", header=TRUE,  ...)
  {
    newname <- paste(dir, name, sep="/")
    read.table(newname, sep=sep, header=header, ...)
  }

########################################################################
## Read in data from GenePix output
## Assuming you include all columns of the output.  
########################################################################

########################################################################/**
# \name{read.genepix}
# 
# \alias{read.genepix}
# 
# \title{Reading in Data Generated by the Image Analysis Package GenePix.}
# 
# \description{
# Reads in a data file in table format and creates a data frame
# with the same number of rows as there are lines in the
#  file, and the same number of columns as there are fields
#  in the file.\cr
#  `read.genepix' reads in the data file generated by the software
#  "GenePix".
# }
# 
# \usage{
# read.genepix(name, dir = ".", sep = "\t", header = T, skip = 26, ...) 
# }
# 
# \arguments{
#  \item{name}{character string naming the data file from which to read the
# data. }
#  \item{dir}{character string naming the directory that contains the
#  data file.} 
# 
#  \item{sep}{the field separator (single  character),  often  "\t"  for
#    tab delimited fields.  If omitted,  any  amount  of white space
#    (blanks or tabs) can separate the fields.  To read  fixed  format
#    files, make  sep  a  numeric vector giving the initial columns of
#    the fields. } 
# 
#  \item{header}{logical flag: if TRUE, then the first line of the
#    file is used as the variable names of the resulting data frame. }
# 
#  \item{skip}{the number of lines of the data file to skip before beginning
#    to read data.}
#  
#  \item{\dots}{parameters for read.table may also be supplied as arguments to
#    the function (see  \code{\link{read.table}}).  }
# }
# 
# \value{as in \code{\link{read.table}}, a data frame
# (\code{\link{data.frame}}) containing a representation of the data
# in the file.  
# }
# 
# \seealso{\code{\link{read.table}}, \code{\link{data.frame}},
# \code{\link{read.spot}}.} 
# 
# \keyword{microarray, GenePix.}
# 
#*/#####################################################################
 
read.genepix <- function (name, dir = ".", sep = "\t", header = TRUE, skip=26, gpname = c("Block", "Col", "Row", "Name", "ID", "X","Y", "Dia", "Rmed", "Rmean", "RSD", "Rbmed", "Rbmean", "RbSD", "Rb1SD", "Rb2SD", "Rbsat","Gmed", "Gmean", "GSD", "Gbmed", "Gbmean", "GbSD", "Gb1SD", "Gb2SD", "Gbsat", "ratiomed", "ratiomean", "medratio", "meanratio", "ratiosd", "Rratio", "RegR2", "FPixels", "BPixels", "summed", "summean", "logratio", "Rmedc", "Gmedc", "Rmeanc", "Gmedc", "flags"),...) 
{
    newname <- paste(dir, name, sep = "/")
    #gpname <- c("Block", "Col", "Row", "Name", "ID", "X","Y", "Dia", "Rmed", "Rmean", "RSD", "Rbmed", "Rbmean", "RbSD", "Rb1SD", "Rb2SD", "Rbsat","Gmed", "Gmean", "GSD", "Gbmed", "Gbmean", "GbSD", "Gb1SD", "Gb2SD", "Gbsat", "ratiomed", "ratiomean", "medratio", "meanratio", "ratiosd", "Rratio", "RegR2", "FPixels", "BPixels", "summed", "summean", "logratio", "Rmedc", "Gmedc", "Rmeanc", "Gmedc", "flags")
    x <- read.table(newname, sep = sep, header = header, skip=skip,quote="", ...)
    colnames(x) <- gpname
    x
}

##########################################################################
# Initialization: slide layout and data matrix of fluorescence intensities
##########################################################################

########################################################################/**
# \name{init.grid}
# 
# \alias{init.grid}
# 
# \title{
# Initialization of Grid Parameters}
# 
# \description{
# Interactive function for specifying the dimensions of the spot
# matrix and the grid matrix. These parameters depend on the printing
# layout of the array, and are used for the within print-tip group
# normalization implemented in \code{\link{stat.ma}} and the spatial
# representation of spot statistics in \code{\link{plot.spatial}}. 
# }
# 
# \usage{
# init.grid()
# }
# 
# \arguments{
# None.
# }
# 
# \value{list containing the following components
# \item{nspot.r}{ the number of rows of spots per grid;}
# \item{nspot.c}{ the number of columns of spots per grid;}
# \item{ngrid.r}{ the number of rows of grids per image;}
# \item{ngrid.c}{ the number of columns of grids per image.}
# }
# 
# \references{
# Spot manual.
# }
# 
# \seealso{
# \code{\link{plot.mva}}, \code{\link{plot.spatial}},
# \code{\link{stat.ma}}, \code{\link{list}}. 
# }
# 
# 
# \author{
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu}
# }
# 
# \examples{
# data(MouseArray)
# # mouse.setup <- init.grid()
# 
# ## Here is what you should see:
# # Enter number of rows of grids per image (ngrid.r): 4
# # Enter number of columns of grids per image (ngrid.c): 4
# # Enter number of rows of spots per grid (nspot.r): 19  
# # Enter number of columns of spots per grid (nspot.c): 21
# # Initialization complete
# }
# 
# \keyword{microarray, grid.}
#*/#########################################################################

init.grid <- function(){
 cat ("Enter number of rows of grids per image (ngrid.r): ")
 ngrid.r <- readline()
 cat ("Enter number of columns of grids per image (ngrid.c): ")
 ngrid.c <- readline()
 cat ("Enter number of rows of spots per grid (nspot.r): ")
 nspot.r <- readline()
 cat ("Enter number of columns of spots per grid (nspot.c): ")
 nspot.c <- readline()
 cat ("Initialization complete\n")
 list(nspot.r = as.integer(nspot.r), nspot.c = as.integer(nspot.c), ngrid.r= as.integer(ngrid.r), ngrid.c = as.integer(ngrid.c))
}          

########################################################################/**
# 
# \name{init.data}
# 
# \alias{init.data}
# 
# \title{Creating a Data Structure for Multi-slide Microarray Experiments}
# 
# \description{
# Interactive function which creates a data structure for multi-slide
# microarray experiments.  The data structure is a list of
# matrices. For each spotted DNA sequence, the list stores raw red and
# green signal intensities as well as red and green background
# intensities. The function also allows the user to add data to an
# existing structure. 
# }
# 
# \usage{
# init.data()
# }
# 
# \arguments{
# None.
# }
# 
# \value{
#   List containing the following components:
#   \item{R}{contains the raw red intensities, R.}
#   \item{G}{contains the raw green intensities, G.}
#   \item{Rb}{contains the background red intensities, Rb.}
#   \item{Gb}{contains the background green intensities, Gb.}
# }
# 
# \author{
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu}
# }
# 
# \examples{
# ## mouse.data <- init.data()
# 
# ## Here is what you should see:
# ## Are you creating a new data matrix or adding new array data
# ## to a prexisting data matrix? 
# ## Enter "n" for creating and "a" for adding new array data: n
# ## Do the names of all your datasets have the following format: 
# ## prefix1, prefix2, prefix3?, ... Here prefix can be any name, 
# ## but the suffixes must be integers 1,2, ..., # of arrays. 
# ## Enter "y" for yes, "n" for no: y
# ## Enter the prefix:mouse
# ## Enter the number of arrays to be processed:6
# ## Enter the name of Cy3 raw data: Gmean
# ## Enter the name of Cy3 background: morphG
# ## Enter the name of Cy5 raw data: Rmean
# ## Enter the name of Cy5 background: morphR
# ## Finished creating new dataset.
# }
# 
# \keyword{microarray.}
# 
#*/#########################################################################

init.data<-function()
{
  ## This file assumes that you have already read in the data.

  cat("Are you creating a new data matrix or adding new array data\n") 
  cat("to a prexisting data matrix? \n")
  cat("Enter \"n\" for creating  and \"a\" for adding new array data: ")
  new.n <- readline()
    
  if(new.n == "a"){
    cat("Enter the name of the existing data matrix: ")
    oname <- readline()
  }
  
  cat("Do the names of all your datasets have the following format: \n")
  cat("prefix1, prefix2, prefix3?, ... Here prefix can be any name, \n")
  cat("but the suffixes must be integers 1,2, ...,  of arrays. \n")
  cat("Enter \"y\" for yes, \"n\" for no: ")
  b.n<-readline()

  if(b.n=="y")
    {
      cat("Enter the prefix:")
      prefixname<-readline()
      cat("Enter the number of arrays to be processed:")
      n<-readline(); n<-as.integer(n)
      dname <- paste(prefixname, 1:n, sep="")
    }
  else if(b.n=="n")
      {
        cat("Enter the number of arrays to be processed:")
        n<-as.integer(readline()); 
        dname<-rep(0,n)

        for(i in 1:n)
          {
            cat(paste("Enter the name of your ", i,"th dataset:"))
            dname[i]<-readline()
          }
      }

  cat ("Enter the name of Cy3 raw data: ")
  name.G <- readline()
  cat ("Enter the name of Cy3 background: ")
  name.Gb <- readline()
  cat ("Enter the name of Cy5 raw data: ")
  name.R <- readline()
  cat ("Enter the name of Cy5 background: ")
  name.Rb <- readline()
  
  if(new.n == "a"){
    res <- eval(as.name(oname))
    action <- "updating"}
  else{
    res <- list(R = NULL, G = NULL, Rb= NULL, Gb=NULL)
    action <- "creating"
  }
  for( i in 1:n){
    tmp <- eval(as.name(dname[i]))[,c(name.R, name.G, name.Rb, name.Gb)]
    res$R <- cbind(res$R, as.numeric(as.vector(tmp[,1])))
    res$G <- cbind(res$G, as.numeric(as.vector(tmp[,2])))
    res$Rb <- cbind(res$Rb, as.numeric(as.vector(tmp[,3])))
    res$Gb <- cbind(res$Gb, as.numeric(as.vector(tmp[,4])))
  }
   
  cat(paste("Finished", action,  "the dataset.\n", sep=" "))
  res
}


########################################################################/**
# 
# \name{init.addinfo}
# 
# \alias{init.addinfo}
# 
# \title{Adding Information to a Data Structure for Multi-slide
# Microarray Experiments} 
# 
# \description{
#   Interactive function which adds other information generated from the
#   output of image analysis software for microarrays to the existing 
#   data structure created using  \code{\link{init.data}}.  
# }
# 
# \usage{
# init.addinfo()
# }
# 
# \arguments{
#   \item{batch}{Character string, this refers to the name of a
#   collection of experiments.} 
#   \item{attri}{Character string, the name of the information to be
#     included in the data structure.  For example, from the output of
#     \tt{Spot}, this argument can be "area", "signal to noise" etc.
#     In other words, these are the column headings from the raw data set.}
#     \item{dataname}{A name of your experimental data.  By default it's
#       named "batch.exp" where batch is the name of the collection of
#       experiments you are interested in.}
#   }
# 
# \value{
#  List containing the following component:
#   
#   \item{R}{contains the raw red intensities, R.}
#   \item{G}{contains the raw green intensities, G.}
#   \item{Rb}{contains the background red intensities, Rb.}
#   \item{Gb}{contains the background green intensities, Gb.}
#   as well as information from other users selected columns information.
# }
# 
# \author{
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
#   Natalie Roberts \email{nroberts@wehi.edu.au}
# }
# 
# \examples{
# ## mouse.data <- init.addinfo("mouse", "area")
# }
# 
# \section{Warning}{The code in the example is not directly executable as
#   it draws upon a particular set of data. This data may be downloaded from
# \url{http://www.stat.berkeley.edu/users/terry/zarray/Software/smacode.html}
# and when loaded appropriately into the user's directory, this example  
# should be executable in its current form. }   
# 
# \keyword{microarray, quality information.}
# 
#*/########################################################################

init.addinfo <- function(batch, attri, dataname=NULL, ...)
{
  if(is.null(dataname)) dataname <- paste(batch, "data", sep=".")
  measure<-NULL 
  nd<-nrow(init.show.exp(batch))
  for(i in 1:nd){
    tmp<-eval(init.read.exp(batch,i, ... ))
    measure<-cbind(measure,tmp[,attri])
  }
  res <- c(eval(as.name(dataname)), list(measure))
  names(res) <- c(names(eval(as.name(dataname))), attri )
  res
}

#########################################################################/**
# 
# \name{init.read.exp}
# \alias{init.read.exp}
# 
# \title{Reads the Output of the Computed Statistics}
# 
# \description{
#   Function displays the 30 measurements computed by the program Spot for
#   each gene in the slide being considered.} 
# 
# \usage{
# init.read.exp(batch, imageid, sep="\t", header=T, ...)}
# 
# \arguments{
#  \item{batch}{batch name of the experiment}
#  \item{imageid}{integer value; the index of the slide which is considered}
#  \item{sep}{the field separator character; the columns of the file
#  will be separated by this character.} 
#  \item{header}{a logical value indicating whether the file contains the
#    names of the variables as its first line.}
#  \item{\dots}{graphical parameters may also be supplied as arguments to
#    the function (see \code{\link{par}}).}
# }
# 
# \value{
#   A matrix containing the 30 columns of computed measurements,
#   corresponding to the rows of different genes in the specified
#   slide. \cr Details regarding these measurements can be found at
#   \url{http://www.cmis.csiro.au/iap/Spot/spotoutput.htm}.}
# 
# \references{Spot manual
# \url{http://www.cmis.csiro.au/iap/Spot/spotmanual.htm}}.
# } 
# 
# \author{Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu}}
# 
# \examples{
# ## apoa1.info <- init.read.exp("apoa1", 1) ## obtains the matrix
# ## of 30 measurements for all the genes spotted on slide 1 of the MouseArray
# ## experiment.}
# 
# \section{Warning}{The code in the example is not directly executable as
#   it draws upon a particular set of data. This data may be downloaded from
#   \url{http://www.stat.berkeley.edu/users/terry/zarray/Software/smacode.html}
# and when loaded appropriately into the user's directory, this
# example should be executable in its current form.} 
# 
# \keyword{measurements, statistics}
#*/#########################################################################

init.read.exp <- function(batch, imageid, sep="\t", header=TRUE, ...)
  {
    tmp<-init.show.exp(batch)
    res<-read.table(tmp[imageid,2], sep=sep, header=header, ...)
    res<-as.matrix(res)
    res<-apply(res, 2, as.numeric)
    res
  }

########################################################################/**
# 
# \name{init.names}
# 
# \alias{init.name.exp}
# \alias{init.show.exp}
# 
# \title{Set and Read the Names of Experimental Data.}
# \description{
#   `init.name.exp' creates a look-up table which contains the names of the
#   experimental data files and the corresponding object names in R. \cr
#   `init.show.exp' displays the look-up table created by
#   \code{init.name.exp}.
# }
# \usage{
# init.name.exp(Robject=F)
# init.show.exp("batch")
# }
# 
# \arguments{
#  \item{Robject}{if TRUE, the function generates a matrix of characters.
#    Otherwise, this matrix is written to a file.}
#  \item{batch}{Character string, this refers to the name of a
#  collection of experiments.}  
# } 
# 
# \value{
#   \code{init.show.exp} returns a list containing the following components:
#   \item{Name in R}{the object names in R;}
#   \item{Filename}{the experimental data filenames, including the full
#     path name for each file.} 
# }
# \references{Spot manual
# \url{http://www.cmis.csiro.au/iap/Spot/spotmanual.htm} 
#   
# \author{Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu}}
# 
# 
# \examples{
# ## init.name.exp() ## To create the look-up table.
# 
# ## This is what you should see:
# ## Are you creating a new batch.exp file or adding new data names
# ## to a prexisting batch.exp file? 
# ## Enter "n" for creating  and "a" for adding new data names: n
# ## Enter the batch name for the new .exp file: mouse1
# ## Enter the number of names of files to be entered: 2
# ##  Enter the R name of your  1 th dataset: m1
# ##  Enter the actual file name including the full path name for m1 ?
# ##  ~/path/image1.data 
# ##  Enter the R name of your  2 th dataset: m2
# ##  Enter the actual file name including the full path name for m2 ?
# ##  ~/path/image2.data 
# ## Finished adding names to .exp file.
# ## NULL
# 
# ## View the look-up table.
# ## init.show.exp("mouse1") 
# ## 
# ##    Name in R           Filename
# ## 1        m1    ~/path/image1.data
# ## 2        m2    ~/path/image2.data
# }
# 
# \keyword{filename}
# 
#*/#########################################################################

init.show.exp <- function(batch)
{
    file <- paste(batch,"exp",sep = ".")
    if (!file.exists(file)){
        stop(paste("File \"", file, "\" does not exist. \n", sep = ""))
	}
    expt <- read.table(file, header = TRUE, as.is = TRUE)
    if (ncol(expt) != 2){
        stop(paste("Should be two columns in experiment name file \"",
            file, "\". \n", sep = ""))
	    }
    colnames(expt) <- c("Name in R", "Filename")
    expt
}

init.name.exp <-function(Robject=FALSE)
{
  ## This file creates the file containing file sources and corresponding 
  ## R names for a batch of experiments
  ## This file assumes that your data exist in the current directory.
  
  cat("\nAre you creating a new batch.exp file or adding new data names\n") 
  cat("to a prexisting batch.exp file? \n")
  cat("Enter \"n\" for creating  and \"a\" for adding new data names: ")
  new.n <- readline()

  while((new.n != "n") & (new.n  != "a")){
     cat("Please enter \"n\" for creating, \"a\" for adding new data names or Ctl-C to quit")
     new.n <- readline()
   }
  
  if(new.n == "a"){
    cat("Enter the batch name of the existing .exp file: ")
    oname <- readline(); 
  }
  
  if(new.n == "n"){
    cat("Enter the batch name for the new .exp file: ")
    oname <- readline(); 
  }
  
  cat("Enter the number of names of files to be entered: ")  
  n<-readline(); n<-as.integer(n)
  dname<-rep(0,n); 
  pname<-rep(0,n);
  
  for(i in 1:n)
    {
      cat(paste("\n Enter the R name of your ", i,"th dataset:"))
      dname[i]<-readline()
      cat(paste("\n Enter the actual file name including the full path name for", dname[i],"?"))
      pname[i]<-readline()
    }

  if(new.n =="n")
    {
      res <- cbind(dname,pname)
      write.table(res, paste(oname, "exp", sep="."),sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
  if(new.n =="a")
    {
      res <-rbind(as.matrix(init.read.exp(oname)), cbind(dname,pname))
      write.table(res, paste(oname, "exp", sep="."),sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
  cat("Finished adding names to .exp file.\n")
  if(!Robject) res <-  NULL
  res
}

########################################################################/**
# \name{init.ctl.index}
# \alias{init.ctl.index}
# \title{Generates co-ordinates of spots.}
# \description{
# Generates the 4 co-ordinates of any spots.
# }
# \usage{
# init.ctl.index(grows, gcols, srows, scols)
# }
# \arguments{
#   \item{grows}{The row index of the grid.}
#   \item{gcols}{The column index of the grid.}
#   \item{srows}{The row index of the spot within the its grid.}
#   \item{scols}{The column index of the spot within the its grid.}
# }
#
# \value{
#   a matrix in which each row contains a vector of 4 integer elements
#   which make up the image coordinates of a gene.
# }
# 
# \examples{
#    x <- init.ctl.index(1:4, 1:4, 1:2, 15:20)
# ## Generates the 4 co-ordinates index to spots in the first 2 rows,
# ## columns 15 to 20 of every print-tips groups.
# }
# \keywords{microarray}
#
# \author{Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu}}
#
# \note{Sorry: No help files yet}
#*/#######################################################################

init.ctl.index <-
function(grows, gcols, srows, scols)
{
ngr <- length(grows)
ngc <- length(gcols)
nsr <- length(srows)
nsc <- length(scols)
t1 <-  rep(grows, rep(nsr * nsc * ngc, ngr) )
t2 <- rep(rep(gcols, rep(nsr * nsc, ngc)), ngr)
t3 <- rep(rep(srows, rep(nsc,nsr)), ngc * ngr)
t4 <- rep(scols, nsr * ngc * ngr)
cbind(t1, t2, t3, t4)
}

##########################################################################
#                                End of file
##########################################################################




































###########################################################################
# Statistics for Microarray Analysis
# Misc. functions
#
# Date : March 19, 2001
#
# Authors: Yee Hwa (Jean) Yang and Jessica Mar
##########################################################################

########################################################################/** 
# \name{is.odd}
# \alias{is.odd}
# \alias{is.even}
# 
# \title{ Determining if a Value is Odd or Even }
# \description{
# A logical flag which determines if a value supplied by the user is
# odd or even. 
# }
# \usage{
# is.odd(x)
# is.even(x)
# }
# 
# \arguments{
#  \item{x}{integer value}
# }
# 
# }
# \value{Logical vectors \code{TRUE} or \code{FALSE} are returned
#   depending on whether the value is odd or even.
# 
# }
# 
# \author{ Jessica Mar }
# 
# \examples{
# is.odd(4)
# ## FALSE
# is.even(100)
# ## TRUE
# }
# 
# \keyword{odd, even}
#*/########################################################################

is.even <- function(x)
  {if(is.numeric(x))
     {if (x %% 2==0) {TRUE}
      else {FALSE}
     }
   
   else{
       print("Warning: Input must be an integer value")
     } 
 }

is.odd <- function(x)
  {if(is.numeric(x))
    {if (x %% 2 == 0) {FALSE}
     else {TRUE}}
 
   else{
       print("Warning: Input must be an integer value")
     }
 }

########################################################################/** 
# \name{id2image}
# 
# \alias{id2image}
# \alias{image2id}
#  
# \title{Converting an id tag to a Set of Image Coordinates and Vice Versa}
# 
# \description{
#  The function \code{id2image} converts an id tag of a gene supplied
#  by a user into a set of image coordinates regarding the location of
#  the gene being 
#  considered. Conversion of image coordinates to an id tag is performed
#  by \code{image2id}.
# }
# 
# \usage{
# id2image(X, layout)
# image2id(x, layout)
# }
# 
# \arguments{
#  \item{X}{an integer value representing the id of a particular gene}
#  \item{layout}{}
#  \item{x}{a vector of 4 integer elements which make up the image
#    coordinates of the gene.}
# }
# \details{
# The image coordinates of a gene correspond to the gene's grid row and
# grid column position within a slide, and the gene's row and column
# position within a grid.  
# }
# 
# \value{
#   \code{id2image} returns a vector of 4 integer elements which is
#   the set of image coordinates. 
#   \code{image2id} returns an integer value which is the gene's id tag.
# }
# 
# \author{Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu}}
# 
# \seealso{\code{\link{MouseArray}}}
# 
# \examples{data(MouseArray)
# # mouse.setup <- init.grid()
# 
# id2image(1024, mouse.setup)
# ## You will see: [1]  1 3 11 16
# ## the grid in which gene 1024 can be found, is in row 1, column 3
# ## and the gene is located in row 11, column 16 of this particular grid.
#  
# image2id(c(2,4,6,8), mouse.setup)
# ## You will see: [1] 2906
# ## the gene located in row 6, column 8 in the grid that is in row 2 and
# ## column 4 is the 2906th gene of the data set.  
# }
# 
# \keyword{image, id}
#*/#######################################################################
 
id2image <- 
function(X, layout)
{
        Grid.row <- layout$ngrid.r; Grid.col <- layout$ngrid.c
        Spot.row <- layout$nspot.r; Spot.col <- layout$nspot.c

        Coord <- rep(0, 4)
        Spot.size <- Spot.row * Spot.col
        ## Calculate Grid row & column coordinates
        Coord[1] <- ((X - 1) %/% (Grid.col * Spot.size)) + 1
        count.Spots <- ((X - 1) %/% Spot.size) + 1
        Coord[2] <- ((count.Spots - 1) %% Grid.col) + 1
        ## Calculate Spot row & column coordinates
        Spot.pt <- X - (count.Spots - 1) * Spot.size
        Coord[3] <- ((Spot.pt - 1) %/% Spot.col) + 1
        Coord[4] <- ((Spot.pt - 1) %% Spot.col) + 1
        Coord
}

image2id <- 
function(x, layout)
{
        Grid.row <- layout$ngrid.r; Grid.col <- layout$ngrid.c
        Spot.row <- layout$nspot.r; Spot.col <- layout$nspot.c

        temp <- Spot.col * Spot.row
        temp * ((x[1] - 1)*Grid.col+(x[2] - 1))+(x[3] - 1)*Spot.col+ x[4]
}



##########################################################################
#                                End of file
##########################################################################
###########################################################################
# Statistics for Microarray Analysis
# MA plots
#
# Date : March 19, 2001
#
# History:
#    May 17, 2001: Incorporate postscript fix from kadowaki@pharmadesign.co.jp
#    March 19, 2001: Some of the plot functions from Rarray.R.
#					
#
# Authors: Sandrine Dudoit and Yee Hwa (Jean) Yang.
##########################################################################

########################################################################/**
# \name{plot.mva}
# 
# \alias{plot.mva}
#   
# \title{M vs. A Plot}
# 
# \description{
# For a single slide, this function produces a scatter plot of log
# intensity ratios \eqn{M = log_2(R/G)} versus average log intensities
# \eqn{A = log_2 \sqrt{RG}}{A = log_2(R*G)/2}, where R and G represent
# the fluorescence intensities in the red and green channels
# respectively. 
# }
# 
# \usage{
# plot.mva(RG, layout, norm="p", pout=T, image.id=1, extra.type="tci",
# crit1=0.025,crit2=crit1, nclass=10, labs=NULL, plot.type="n",
# col.ex=NULL, ...) 
# }
# 
# \arguments{
#   \item{RG}{A list with at least 4 elements.  Each element of the list
#     being a matrix with p rows for p genes and n columns for n slides. 
#     The first element 'R' contains the raw red intensities,
#     the second element 'G' contains the raw green intensities,
#     the third element 'Rb' contains the background red intensities and
#     the 4th element 'Gb' contains the background green intensities.
#     This data structure can be generated by an interactive function
#     \code{\link{init.data}}.}  
# 
#   \item{layout}{a list specifying the dimensions of the spot matrix
#   and the grid matrix.  This can be generated by calling
#   \code{\link{init.grid}}.} 
# 
#   \item{norm}{character string, one of "n", "m", "l", "p" or "s".  This
#     argument  specifies the type of normalization method to be
#     performed: "n" no normalization between the 2 channels; "m"
#     \code{\link{median}} normalization, which sets the median of log
#     intensity ratios to zero; "l" global \code{\link{lowess}}
#     normalization; "p" print-tip group lowess normalization and "s"
#     scaled print-tip group lowess normalization.} 
#   
#   \item{pout}{if TRUE, an M vs. A plot will be produced.  Otherwise, the
#     function returns the normalized log intensity ratios M and the mean
#     log intensities A for each gene.} 
#   
#   \item{image.id}{integer value; the index of the slide which is considered.}
# 
#   \item{extra.type}{a character string, one of "t", "p", "tci","pci" or
#     "lci".  This argument specifies the type of plot to be drawn.  The
#     possible types are: \cr
#     * "t" for text, \cr
#     * "p" for points, \cr
#     * "tci" for highlighting a certain proportion of extreme `M' values
#     by text,\cr
#     * "pci" for highlighting a certain proportion of extreme `M' values
#     by points,\cr
#     * "lci" for including 2 intensity dependent lines where a 
#     prespecified proportion of points have more extreme `M' values. 
#   }
#   \item{crit1}{The number of points to be highlighted on the M vs A
#     plot.  If crit1 < 1, the crit1*100\% spots with the smallest M
#     values will be highlighted. If crit1 >= 1, the crit spots 
#    with the smallest M values are highlighted.} 
#  \item{crit2}{Similar to "crit1".   If crit2 < 1, the crit2*100\% spots
#    with the largest M values will be highlighted. If crit2 >= 1, the
#    crit2 spots with the smallest M values are highlighted.}
#  \item{nclass}{A single number giving the approximate number of
#    intensity dependent groups to consider.}
#  \item{labs}{one or more character strings or expressions specifying the
#    text to be written.  If this string is not specified, by
#    default the index of the vector `M' will be used.}
#  \item{plot.type}{a character string, this argument is either "n", "r"
#    or "b".  The different number of plots to be included are:\cr
#    * "n" for normalised M vs A plot, \cr
#    * "r" for unnormalised M vs A plot, and \cr 
#    * "b" both unnormalised and normalised M vs A plots.
#   }
#  \item{col.ex}{The colour used for the highlighting extreme points,
#  lines or text.} 
#  \item{\dots}{graphical parameters may also be supplied as arguments to the
#     function (see \code{\link{par}}).  }
# }
# 
# \value{A plot is created on the current graphics device.  The top
# plot is based on unnormalized log ratios and the bottom plot is
# based on normalized log ratios.} 
# 
# \details{M vs. A plots tend to be more revealing than their log R
# vs. log G counterparts in terms of identifying spot artifacts and
# detecting intensity dependent patterns in the log ratios. They are
# also very useful for normalization.} 
# 
# \references{S. Dudoit, Y. H. Yang, M. J. Callow, and T. P. Speed. Statistical
# methods for identifying differentially expressed genes in replicated
# cDNA microarray experiments (Statistics, UC Berkeley, Tech Report \#
# 578).}
# 
# \author{
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} \cr
#   Natalie Roberts, \email{nroberts@wehi.edu.au}
# }
# 
# \seealso{\code{\link{ma.func}}, \code{\link{plot.smooth.line}},
# \code{\link{stat.ma}}, \code{\link{lowess}}, \code{\link{plot}}.} 
# 
# \examples{
# data(MouseArray)
# # mouse.setup <- init.grid()
# # mouse.data <- init.data() ## see \emph{init.data}
# mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# 
# ## Look at the normalized second data sets in the list using points to
# ## highlight large positive or large negative ratios.
# plot.mva(mouse.data, mouse.setup, norm="l", 2, extra.type="pci",
# plot.type="n") 
# 
# ## Look at the both unnormalized and normalized first data sets in the
# ## list using text to highlight large positive or negative ratios.
# ## plot.mva(mouse.data, mouse.setup, norm="l", 2, extra.type="tci", plot.type="b") 
# }     
# 
# \keyword{microarray}
# 
#*/########################################################################

plot.mva <- function(x, layout, norm="p", pout=TRUE, image.id=1, extra.type="tci", crit1=0.025,crit2=crit1, nclass=10, labs=NULL, plot.type="n", col.ex=NULL, pch=".", ...)
{
#  RG <- x
  ma.func(R = x$R[,image.id], G=x$G[,image.id], Rb=x$Rb[,image.id], Gb = x$Gb[,image.id], layout=layout, norm=norm, pout=pout, extra.type=extra.type, crit1=crit1, crit2=crit2, nclass=nclass, labs=labs, plot.type=plot.type, col.ex=col.ex, pch=pch,...)
}

########################################################################/**
# \name{plot.smooth.line}
# 
# \alias{plot.smooth.line}
# 
# \title{Adding Lowess Lines to Current Plot}
# 
# \description{
#  This function adds a \code{\link{lowess}} line to the current
#  plot.  The  type of line can be specified as well as other
#  parameters.} 
# 
# \usage{
# plot.smooth.line(A, M, f=0.1, ...)
# }
# 
# \arguments{
#  \item{A}{a vector giving the x-coordinates of the points in the scatter
#           plot. In the microarray context, this could be a vector of
#           average log intensities.} 
# 
#  \item{M}{a vector giving the y-coordinates of the points in the scatter
#           plot. In the microarray context, this could be a vector of
#           log intensity ratios.} 
# 
#  \item{f}{the smoother span. This gives the proportion of points in the
#           plot which influence the smoothness at each value. Larger
# 	  values give greater smoothness. } 
# 
# \item{\dots}{graphical parameters may also be supplied as arguments
# to the function (see \code{\link{par}}).}        
# }
# 
# \value{
#  lines are added to the current plot.
# }
# 
# \note{An M vs A plot must be constructed \bold{prior} to the execution of this function.}
# 
# \references{ Chambers, J. M., Cleveland, W. S., Kleiner, B. and Tukey,
# P. A. (1983). Graphical Methods for Data Analysis. Wadsworth, Belmont,
# California. }
# 
# \seealso{ \code{\link{plot.mva}}, \code{\link{stat.ma}},
# \code{\link{lines}}, \code{\link{lowess}}, \code{\link{smooth}}. 
# }
# 
# \examples{
# data(MouseArray)
# ## mouse.setup <- init.grid()
# ## mouse.data <- init.data()
# 
# plot.mva(mouse.data, mouse.setup)
# plot.smooth.line(mouse.lratio$A, mouse.lratio$M)
# }
# 
# \keyword{microarray, lowess.}
#*/########################################################################

plot.smooth.line  <- function(x, M, f = 0.1, ...)
{
#  A <- x
  ind <- !(is.na(x) | is.na(M) | is.infinite(x) | is.infinite(M))
  #lines(lowess(A[ind], M[ind], f = f), ...)
  lines(approx(lowess(x[ind], M[ind], f = f)), ...)  
}

########################################################################/**
# \name{plot.confband.lines}
# 
# \alias{plot.confband.lines}
# 
# \title{Adding Lines Satisfying a Confidence Criterion to the Current M
#   vs A Plot}
# 
# \description{
#   This function adds 2 lines outlining the pointwise (intensity
#   dependent) confidence band on the M vs A plot.  The lines are drawn
#   such that a prespecified proportion of points are outside the 2
#   confidence curves.
#   The type of line may be specified as well as other parameters.}
# 
# \usage{
# plot.confband.line(A, M, crit1=0.025, crit2=crit1, nclass=10, ...)
# }
# 
# \arguments{
#  \item{A}{a vector giving the x-coordinates of the points in the scatter
#           plot. In the microarray context, this could be a vector of
#           average log intensities.} 
# 
#  \item{M}{a vector giving the y-coordinates of the points in the scatter
#    plot. In the microarray context, this could be a vector of log
#    intensity ratios.} 
# 	
#  \item{crit1}{The proportion of points less than the lower confidence
#    curve.  This takes a decimal value between 0 and 1. }
#  \item{crit2}{The proportion of points greater than the upper confidence
#    curve.  By default, this has the same value as "crit1".}
#  \item{nclass}{A single number giving the approximate number of
#  intensity depedent groups to consider.} 
#  \item{\dots}{graphical parameters may also be supplied as arguments
#  to the function (see  \code{\link{par}}).}        
# }
# 
# \value{
#   Lines are added to the current plot.
# }
# 
# \note{
#   An M vs A plot must be constructed \bold{prior} to the execution
#   of this function.} 
# 
# \seealso{ \code{\link{plot.mva}}, \code{\link{stat.ma}},
#   \code{\link{lines}}, \code{\link{matlines}},
#   \code{\link{plot.confband.text}}, \code{\link{plot.confband.points}} .
# }
# 
# \examples{data(MouseArray)
# ## mouse.setup <- init.grid
# ## mouse.data <- init.data
# 
# ## To display an M vs A plot of the data 
# plot.mva(mouse.data, mouse.setup) 
# 
# ## Calculate M and A values 
# mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# 
# ## To add default upper and lower confidence curves line to the M vs A plot
# plot.confband.lines(mouse.lratio$A, mouse.lratio$M) 
# }
# 
# \keyword{microarray, point-wise confidence band.}
#*/########################################################################

plot.confband.lines<-function (x, M, crit1=0.025, crit2=crit1, nclass=10, ...)
{
#  A <- x
  if (crit1 >= 1) crit1 <- crit1 / length.na(M)
  if (crit2 >= 1) crit2 <- crit2 / length.na(M)
  cutoff<-NULL
  Abin <- quantile(x, probs=seq(0, nclass, 1)/nclass, na.rm=TRUE)
  for(i in (1:nclass) ){
    tmpind<-(Abin[i]<=x)&(x<Abin[i+1])
    xtmp <- M
    xtmp[!tmpind]<-NA
    n1<-sum.na(tmpind)
    cutoff <- rbind(cutoff,quantile.na(xtmp, probs=c(crit1, (1-crit2))))
  }
  matlines(Abin[-1],cutoff, ... )
}

########################################################################/**
# \name{plot.confband.points}
# 
# \alias{plot.confband.points}
# 
# \title{Highlights a Set of Points on the Current M vs A Plot}
# 
# \description{
#   This function highlights a prespecified proportion of extreme points
#   on the M vs A plots. 
# }
# 
# \usage{
# plot.confband.points(A, M, crit1=0.025, crit2=crit1, nclass=10, ...)
# }
# 
# \arguments{
#  \item{A}{a vector giving the x-coordinates of the points in the scatter
#           plot. In the microarray context, this could be a vector of
#           average log intensities.} 
# 
#  \item{M}{a vector giving the y-coordinates of the points in the scatter
#    plot. In the microarray context, this could be a vector of log
#    intensity ratios.} 
# 	
#  \item{crit1}{The number of points to be highlighted on the M vs A plot.
#    If crit1 < 1, the crit1*100\% spots with the smallest M values
#    will be highlighted. If crit1 >= 1, the crit spots with the
#    smallest M values are highlighted.}  
#  \item{crit2}{Similar to "crit1".   If crit2 < 1, the crit2*100\%
#  spots with the largest M values will be highlighted. If crit2 >= 1,
#  the crit2 spots with the smallest M values are highlighted.}  
#  \item{nclass}{A single number giving the approximate number of
#  intensity depedent groups to consider.} 
#  \item{\dots}{graphical parameters may also be supplied as arguments
#  to the function (see \code{\link{par}}).}        
# }
# 
# \value{
#   Points are added to the current plot.
# }
# 
# 
# \seealso{\code{\link{plot.mva}}, \code{\link{stat.ma}},
#   \code{\link{lines}}, \code{\link{matlines}},
#   \code{\link{plot.confband.text}}, \code{\link{plot.confband.lines}} .
# }
# 
# \note{An M vs A plot must be constructed \bold{prior} to the
# execution of this function.} 
# 
# \examples{data(MouseArray)
# ## mouse.setup <- init.grid()
# ## mouse.data <- init.data()
# 
# plot.mva(mouse.data, mouse.setup) ## an M vs A plot 
# 
# mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# 
# plot.confband.points(mouse.lratio$A, mouse.lratio$M)
# 
# ## 2.5\% of the spots with the smallest and largest M values are
# ## highlighted on the M vs A plot. 
# }
# 
# \keyword{microarray, point-wise confidence band.}
# 
#*/########################################################################
 

plot.confband.points<-function (x, M, crit1=0.025, crit2=crit1, nclass=10, col.ex=NULL, ...)
{
  ## quantile.na removes infinite too...quantile(x, na.rm=F) doesn't.

  # A <- x
  
  if (crit1 >= 1) crit1 <- crit1 / length.na(M)
  if (crit2 >= 1) crit2 <- crit2 / length.na(M)  
  txtA<-(rep(FALSE,length(x)))
  Abin <- quantile(x, probs=seq(0, nclass, 1)/nclass, na.rm=TRUE)
  for(i in 1:nclass){
    tmpind<-(Abin[i]<=x)&(x<Abin[i+1])
    xtmp <- M
    xtmp[!tmpind]<-NA
    n1<-sum.na(tmpind)
    cutoff <- quantile.na(xtmp, probs=c(crit1, (1-crit2)))
    vals<- ((xtmp < cutoff[1]) | (xtmp > cutoff[2]))
    txtA[vals]<-TRUE 
  }
  points(x[txtA],M[txtA],pch=18, col=col.ex,...)
}


########################################################################/**
# \name{plot.confband.text}
# 
# \alias{plot.confband.text}
# 
# \title{Add Selected Text to an M vs A Plot}
#   
# \description{`text' draws the strings given in the vector `labs' at the
#   coordinates given by `M' and `A'}
# 
# \usage{
# plot.confband.text(A, M, crit1=0.025, crit2=crit1, nclass=10,
# labs=NULL, output=F, ...)  
# }
# 
# \arguments{
#  \item{A}{a vector giving the x-coordinates of the points in the scatter
#           plot. In the microarray context, this could be a vector of
#           average log intensities.} 
#  
#  \item{M}{a vector giving the y-coordinates of the points in the scatter
#    plot. In the microarray context, this could be a vector of log
#    intensity ratios.} 
# 	
#  \item{crit1}{The number of points to be highlighted on the M vs A plot.
#    If crit1 < 1, the crit1*100\% spots with the smallest M values
#    will be highlighted. If crit1 >= 1, the crit spots  with the
#    smallest M values are highlighted.}  
#  \item{crit2}{Similar to "crit1".   If crit2 < 1, the crit2*100\%
#  spots with the largest M values will be highlighted. If crit2 >= 1,
#  the crit2 spots with the largest M values are highlighted.}  
#  \item{nclass}{A single number giving the approximate number of
#  intensity depedent groups to consider.} 
#  \item{labs}{ one or more character strings or expressions specifying the
#    text to be written.  If this string is not specified, by
#    default the index of the vector `M' will be used.}
#  \item{output}{logical, defaulting to `FALSE'. If `TRUE' a vector
#    containning the index to the vector `M' that are  
#    highlighted.} 
#  \item{\dots}{graphical parameters may also be supplied as arguments
#  to the function (see \code{\link{par}}).}        
# }
# 
# \note{An M vs A plot must be constructed \bold{prior} to the execution of this function.}
# 
# \examples{data(MouseArray)
# ## mouse.setup <- init.grid()
# ## mouse.data <- init.data()
# 
# plot.mva(mouse.data, mouse.setup) ## an M vs A plot
# 
# mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# 
# plot.confband.text(mouse.lratio$A, mouse.lratio$M)
# ## 2.5\% of the spots with the largest and smallest M values are
# ## highlighted on the M vs A plot, and each spot is assigned the
# ## default label of its corresponding index value. 
# }
# 
# \seealso{ \code{\link{plot.mva}}, \code{\link{stat.ma}},
#   \code{\link{lines}}, \code{\link{matlines}},
#   \code{\link{plot.confband.lines}}, \code{\link{plot.confband.points}} .
# }
# 
# \keyword{microarray, point-wise confidence band.}
#*/########################################################################
 

 plot.confband.text <- 
function (x, M, crit1=0.025, crit2=crit1, nclass=10, labs=NULL, output=FALSE, ...) 
{
#  A <- x
  if (crit1 >= 1) crit1 <- crit1 / length.na(M)
  if (crit2 >= 1) crit2 <- crit2 / length.na(M)

  txtA<-(rep(FALSE,length(x)))
  Abin <- quantile.na(x, probs=seq(0, nclass, 1)/nclass)
  for(i in 1:nclass){
    tmpind<-(Abin[i]<=x)&(x<Abin[i+1])
    xtmp <- M    
    xtmp[!tmpind]<-NA
    n1<-sum.na(tmpind)
    cutoff <- quantile.na(xtmp, probs=c(crit1, (1-crit2)))
    vals<- ((xtmp < cutoff[1]) | (xtmp > cutoff[2]))
    txtA[vals]<-TRUE
  }
  if(is.null(labs)) labs <- as.character(1:length(M))
  text(x[txtA],M[txtA],labels=labs[txtA], ...)
  if(output)
    res <- txtA
  else res <- NULL
  res
}


##########################################################################
#                                End of file
##########################################################################



###########################################################################
# Statistics for Microarray Analysis
# Function dealing with NA's
#
# Date : March 19, 2001
#
# Authors: Sandrine Dudoit and Yee Hwa (Jean) Yang.
##########################################################################

##########################################################################
# Basic statistics functions that are able to handle missing values
##########################################################################

########################################################################/**
# \name{na}
# 
# \alias{log.na}
# \alias{sum.na}
# \alias{mean.na}
# \alias{var.na}
# \alias{cor.na}
# \alias{quantile.na}
# \alias{length.na}
# \alias{order.na}
# \alias{scale.na}
# \alias{prod.na}
# 
# \title{Basic Statistical Functions for Handling Missing Values}
# 
# \description{
# Basic statistical functions for handling missing values or NA. \cr 
# In \code{log.na}, \code{sum.na}, \code{mean.na} and \code{var.na},
# \code{quantile.na}, \code{length.na}, missing values are omitted
# from the calculation. \cr 
# The function \code{cor.na} calls \code{cor} with the argument
# \code{use="pairwise.complete.obs"}. \cr 
# The function \code{order.na} only handles vector arguments and not
# lists.  However, it gives the option of omitting the NAs
# (\code{na.last=NA}), of placing the NAs at the start of the ordered
# vector (\code{na.last=F}) or at the end (\code{na.last=T}). \cr 
# The function \code{scale.na} is a modified version of
# \code{\link{scale}} which allows NAs in the variance calculation. If
# \code{scale = T}, the function \code{f} in \code{scale.na} uses
# \code{var.na} to perform the variance calculation.
# The function \code{prod.na} is similar to the \code{\link{prod}}
# function with \code{na.rm=TRUE}. This function returns the product of
# all the values present in its arguments, omitting any missing values.
# }
# 
# \author{
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu}
# }
# 
# \seealso{\code{\link{log}}, \code{\link{sum}}, \code{\link{mean}},
#   \code{\link{var}}, \code{\link{cor}}, \code{\link{order}},
#   \code{\link{scale}}, \code{link{prod}}.}
# 
# \keyword{log, sum, mean, variance, correlation, order, scale,
# product, missing values, NA.} 
# 
#*/#########################################################################

 
mean.na <- function(x,...)
{
        mean(x[!(is.na(x) | is.infinite(x))])
}
 
 var.na <- function(x)
{
        res <- NA
        tmp <- !(is.na(x) | is.infinite(x))
        if(sum(tmp) > 1)
                res <- var(x[tmp])
        res
}

cor.na <- function(x)
{
  cor(x, use="pairwise.complete.obs")
}

sum.na <- function(x)
{
        res <- NA
        tmp <- !(is.na(x) | is.infinite(x))
        if(sum(tmp) > 0)
                res <- sum(x[tmp])
        res
}


length.na <- function(x, ...)
{
   tmp <- !(is.na(x) | is.infinite(x))
   length(x[tmp],...)
 }

log.na <- function(x, ...)
{
  log(ifelse(x > 0, x, NA), ...)
}


quantile.na <- function(x, ...)
 {          
   tmp <- !(is.na(x) | is.infinite(x))
   quantile(x[tmp],...)
 }
   
order.na <- function (x, na.last = TRUE) 
{
    y <- order(x)
    n <- sum(is.na(x))
    tmp <- (length(x) - n + 1):length(x)
    if (!is.na(na.last)) {
        if (na.last) 
            res <- y
        if (!na.last)
          {
            if(n == 0)
              res <- y
            else
              res <- c(y[tmp], y[-tmp])
          }
      }
    if (is.na(na.last)) {
        warning("NA's discarded")
        res <- y[-tmp]
    }
    res
}

scale.na<-function(x, center = TRUE, scale = TRUE)
{
  x <- as.matrix(x)
  nc <- ncol(x)

  if (is.logical(center)) {
     if (center)
       x <- sweep(x, 2, applyy(x, 2, mean, na.rm=TRUE))
    }
  else if (is.numeric(center) && (length(center) == nc))
    x <- sweep(x, 2, center)
  else
    stop("Length of center must equal the number of columns of x")
  
  if (is.logical(scale)) {
    if (scale) {
      f <- function(v) {
        sqrt(var.na(v))
      }
      x <- sweep(x, 2, apply(x, 2, f), "/")                   
    }
    }
  else if (is.numeric(scale) && length(scale) == nc)
    x <- sweep(x, 2, scale, "/")
  else
    stop("Length of scale must equal the number of columns of x")
    x
}

prod.na <- function (x) 
{
  prod(x[!(is.na(x) | is.infinite(x))])
}


##########################################################################
#                                End of file
##########################################################################
###########################################################################
# Statistics for Microarray Analysis
# Misc plots
#
# Date : March 19, 2001
#
# History:
#    March 19, 2001: Some of the plot functions from Rarray.R.
#					
#
# Authors: Sandrine Dudoit and Yee Hwa (Jean) Yang.
##########################################################################

##########################################################################
# Exploratory plots for single slide
##########################################################################

########################################################################/**
# \name{plot.svb}
# 
# \alias{plot.svb}
# \alias{svb.func}
# 
# \title{Plot of Signal vs. Background}
# 
# \description{
# Produces a scatter plot of background corrected signal intensities
# and background intensities.
# }
# 
# \usage{
# plot.svb(x, channel="R", image.id=1, S.isbgcor=F, ...)
# }
# 
# \arguments{
#   \item{x}{a numeric list of signal and background intensities, can
#   be raw or background corrected data.}   
# 
#  \item{channel}{the specific channel to which the intensities to be
#    considered, correspond to, that is, either red or green. The default
#    channel is red.}
# 
#  \item{image.id}{integer value; the index of the slide which is considered}
# 
#  \item{S.isbgcor}{logical flag, equal to TRUE if the signal intensities in
#    x contain background corrected signal intensities instead of raw
#    signal intensities. By default this is set to FALSE.}
# 
#  \item{\dots}{graphical parameters may also be supplied as arguments
#  to the function (see \code{\link{par}}).} 
# 
# }
# 
# \value{a plot is created on the current graphics device.}
# 
# \author{
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} \cr
#   Jessica Mar}
# }
#   \seealso{\code{\link{plot}}.} 
# 
# \examples{
# data(MouseArray)
# # mouse.setup <- init.grid()
# # mouse.data <- init.data() ## see \emph{init.data}
# 
# plot.svb(mouse.data, "green", 3) 
# ## thiscreates a plot of the signal versus background intensities 
# ## for the green channel, using data collected from the third slide. 
# }
# 
# \keyword{microarray, background}
# 
#*/######################################################################/**
plot.svb <- function(x, channel="R", image.id=1, S.isbgcor=FALSE, ...)
   {if(is.list(x))
     if ((channel=="R") | (channel=="r") | (channel=="red"))
     {svb.func(x$R[, image.id], x$Rb[, image.id], S.isbgcor=FALSE, ...)}
       
     if ((channel=="G") | (channel=="g") | (channel=="green"))
     {svb.func(x$G[, image.id], x$Gb[, image.id], S.isbgcor=FALSE, ...)}
   }

svb.func <- function(Signal, Bg, S.isbgcor = FALSE, ...)
{
  if(S.isbgcor){
    ind <- log.na(Signal, 2) < quantile(log.na(Signal, 2), 0.75,  na.rm=TRUE)
    plot(log.na(Bg,2)[ind], log.na(Signal,2)[ind], 
	 xlab="Background", ylab="Signal", ...)
  }
  if(!S.isbgcor){
    ind <- log.na(Signal-Bg, 2) < quantile(log.na(Signal-Bg, 2), 0.75,  na.rm=TRUE)
    plot(log.na(Bg,2)[ind], log.na(Signal-Bg,2)[ind], 
	 xlab="Background", ylab="Signal", ...)
  }
}

#######################################################
# plot.print.tip.lowess - a function to print individual
# lowess curves for each print tip super imposed onto an
# m vs a  plot.
#
# TO be added a linetype palette, functionality to add
# an index
#
#######################################################

plot.print.tip.lowess <- function (x, layout, norm = "n", image.id = 1, palette = rainbow(layout$ngrid.r*layout$ngrid.c), lty.palette = rep(1,layout$ngrid.r*layout$ngrid.c), ...) 
{
    tmp <- ma.func(R = x$R[, image.id], G = x$G[, image.id], 
        Rb = x$Rb[, image.id], Gb = x$Gb[, image.id], layout, 
        norm = norm, pout = FALSE, ...)
    plot(tmp$A, tmp$M, xlab = "A", ylab = "M", ...)
    npin <- layout$ngrid.r * layout$ngrid.c
    nspot <- layout$nspot.c * layout$nspot.r
    pin <- rep(1:npin, rep(nspot, npin))
    for (i in 1:npin) {
        index <- pin == i
        tM <- tmp$M[index]
        tA <- tmp$A[index]
        ind2 <- is.na(tM) | is.na(tA) | is.infinite(tM) | is.infinite(tA)
        smoothnum <- lowess(tA[!ind2], tM[!ind2])
        lines(approx(smoothnum), col = palette[i], lty = lty.palette[i], ...)
    }
}





##########################################################################
# Diagnostic plots for multiple slides only
##########################################################################

########################################################################/**
# \name{plot.qq}
# 
# \alias{plot.qq}
# 
# \title{ Histogram and Normal Quantile-Quantile plot}
# 
# \description{Produces a histogram and a normal Quantile-Quantile plot
# of the data. The points corresponding to genes with statistics
# less/greater than a user defined threshold are highlighted. The
# histogram and Q-Q plots are displayed on the same page. 
# }
# 
# \usage{
# plot.qq(x, name, low=-5, high=5)
# }
# 
# \arguments{
# 
#  \item{x}{a numeric vector containing the statistics whose histogram
#  and Q-Q plot will be produced. Missing values (NAs) are allowed.}
# 
#  \item{name}{title for the plots.}
# 
#  \item{low}{lower threshold: points with statistic < low are colored
#  in green.}
# 
#  \item{high}{upper threshold: points with statistic > high are
#  colored in red.} 
# } 
# 
# \references{Chambers, J. M., Cleveland, W. S., Kleiner, B. and
# Tukey, P. A. (1983). Graphical Methods for Data Analysis. Wadsworth,
# Belmont, California. 
#  
# Hoaglin, D. C., Mosteller, F. and Tukey, J.  W., editors
#  (1983). Understanding Robust and Exploratory Data Analysis. Wiley,
#  New York.         
# }
# 
# \author{
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} \cr
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} 
# }
# 
# \seealso{\code{\link{plot.spatial}}, \code{\link{plot.t2}},
# \code{\link{stat.t2}}, \code{\link{hist}}, \code{\link{qqnorm}}.} 
# 
# \examples{
# data(MouseArray)
# ## mouse.setup <- init.grid()
# ## mouse.data <- init.data() ## see \emph{init.data}
# ## mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# 
# ## Calculation of t-statistics
# ## cl <- c(rep(1,3), rep(2,3))
# ## mouse.t2 <- stat.t2(mouse.lratio, cl)
# 
# ## Diagnostic plots
# plot.qq(mouse.t2$t, "Mouse")
# }     
# 
# \keyword{microarray, histogram, qqplot.}
# 
#*/######################################################################/**

plot.qq <- function(x,name, low=-5, high=5,...)
{
  par(mfrow=c(2,1))
  hist(x,xlab="t",nclass=100,main=paste(name," Histogram and quantile-quantile plot of t-statistics", sep=":"),col=9,cex=0.8)
  tmp<-qqnorm(x,plot=FALSE)
  plot(tmp,pch=".",xlab="Quantiles of standard normal",ylab="t")
  points(tmp$x[tmp$y<low],tmp$y[tmp$y<low],pch="*",col=6,cex=1.2)
  points(tmp$x[tmp$y>high],tmp$y[tmp$y>high],pch="*",col=2,cex=1.2)
}

 
#########################################################################/**
# \name{plot.qqline}
# 
# \alias{plot.qqline}
# 
# \title{Add Line Going Through the Quantiles of a Q-Q Plot}
# 
# \description{
# This function adds a line to a quantile-quantile plot which passes
# through user defined quantiles. This function is similar to, but
# more general than, \code{\link{qqline}} because the reference
# distribution need not be the standard normal distribution and the
# quantiles need not be the first and third quartiles. \cr 
# Graphical parameters may be given as arguments to \code{plot.qqline}. 
# }
# 
# \usage{
# plot.qqline(x, y, a=0.25, ...)
# }
# 
# \arguments{
# \item{x}{the reference (first) sample for the Q-Q plot, for a normal Q-Q
#   plot this would be the quantiles of a N(0,1) random sample.}
# \item{y}{the data.}
# }
# \item{a}{a number between 0 and 1. A line is drawn which connects the
#   \code{a} and \code{1-a} quantile points. The default line passes
#   through the first and third quantiles.}
#  \item{\dots}{graphical parameters may also be supplied as arguments
#  to the function (see \code{\link{par}}).} 
# }
# 
# \seealso{\code{\link{qqplot}}, \code{\link{qqnorm}}, \code{\link{qqline}}.  }
# 
# \examples{
# data(MouseArray)
# # mouse.setup <- init.grid()
# # mouse.data <- init.data() ## see \emph{init.data}
# # mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# 
# ## Calculation of t-statistics
# ## cl <- c(rep(1,3), rep(2,3))
# ## mouse.t2 <- stat.t2(mouse.lratio, cl)
# 
# ## Diagnostic plots
# plot.qq(mouse.t2$t, "Mouse")
# 
# ## Using the QQline function
# q <- quantile(rnorm(1000))
# plot.qqline(q, mouse.t2$t)
# }    
# 
# \keyword{Q-Q plots, quartiles}
# 
#*/#########################################################################

plot.qqline<-function(x,y,a=0.25, ...)
{
    y <- quantile(y[!is.na(y)],c(a, 1-a))
    x <- quantile(x[!is.na(x)],c(a, 1-a))
    points(x,y,...)
    slope <- diff(y)/diff(x)
    int <- y[1]-slope*x[1]
    abline(int, slope, ...)
}                       

#########################################################################/**
# \name{plot.scale.box}
# \alias{plot.scale.box}
# 
# \title{Box plots for microarray}
# \description{
# Produce box-and-whisker plot(s) of the given (grouped) values.
# }
# \usage{
# plot.scale.box(x, layout, x.names=NULL, ...)
# }
# 
# \arguments{
#   \item{x}{a vector or a matrix.}
#   \item{layout}{ a list specifying the dimensions of the spot matrix
#    and the grid matrix.  This can be generated by calling
#    \code{\link{init.grid}}.}
#   \item{x.names}{group labels which will be printed under each boxplot.}
#   \item{\dots}{further arguments to the default boxplot method and graphical
#     parameters may also be passed as arguments, see \code{\link{par}}.}
# }
# \details{
#   If x is a vector, this function will produce n boxplots where n is
#   number of print-tips groups.   If x is a matrix, this function will
#   produce n boxplots where n is number of columns in the matrix.  
# }
# 
# \author{Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu}}
# 
# \seealso{\code{\link{boxplot}}, \code{\link{bxp}}}
# 
# \examples{
#      data(MouseArray)
#      # mouse.setup <- init.grid() 
#      # mouse.data <- init.data() ## see \emph{init.data} 
#      mouse.lratio <- stat.ma(mouse.data, mouse.setup)
#      ## Producing boxplots for different print-tips groups.
#      plot.scale.box(mouse.lratio$M[,1], mouse.setup)
#      ## Producing boxplots for different slides.
#      plot.scale.box(mouse.lratio$M)
# }
# \keyword{boxplots, microarray}
# 
#*/#########################################################################


plot.scale.box <- 
function(x, layout=NULL, x.names=NULL, ...)
  {
    n <- layout$nspot.r * layout$nspot.c * layout$ngrid.r * layout$ngrid.c
    nperpin <- layout$nspot.r * layout$nspot.c
    npin <- layout$ngrid.r * layout$ngrid.c
    
    if(is.vector(x)){
      if((length(x) != n) & (is.null(x.names)))
        {
          stop(" Error: Length of vector different from total number of spots and vector has no row.name.\n")
        }
      if ((length(as.vector(x)) != n) & (!is.null(x.names)))
        {
          y <- x; x <- rep(NA, n);
          x[as.integer(x.names)] <- y
      }
      xmat <- matrix(x, nrow = nperpin)
       vect <- TRUE
    }

    if(is.matrix(x))
      xmat <- x
    
    boxplot(data.frame(xmat), names=x.names, ...)
  }



########################################################################/**
# \name{plot.t2}
# 
# \alias{plot.t2}
# 
# \title{Diagnostic Plots for Two-Sample t-statistics}
# 
# \description{
# Plots of two-sample t-statistics, |t-numerator| and
# t-denominator against average A, and plot of |t-numerator| against
# t-denominator. For each spot on a given slide, \eqn{A = log_2
# \sqrt{RG}}{A = log_2(R*G)/2}, where (R,G) denotes the red and green
# fluorescence intensity pair. Points with t-statistics exceeding user
# defined thresholds are highlighted. 
# } 
# 
# \usage{
# plot.t2(X, main.title="T plots", low=-5, high=5)
# }
# 
# \arguments{
# 
#  \item{X}{output from the function \code{\link{stat.t2}}.}
# 
#  \item{main.title}{title for the plot.}
# 
#  \item{low}{lower threshold for t-statistic: points with t<low are
#  colored in green.} 
# 
#  \item{high}{upper threshold for t-statistic: points with t>high are
#  colored in red.} 
# }
# 
# \author{
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} \cr
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} 
# }
# \seealso{\code{\link{stat.t2}}, \code{\link{t2stat.func}},
# \code{\link{plot}}, \code{\link{t.test}}.} 
# 
# \examples{
# data(MouseArray)
# # mouse.setup <- init.grid()
# # mouse.data <- init.data() ## see \emph{init.data}
# # mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# 
# ## Calculation of t-statistics
# ## cl <- c(rep(1,3), rep(2,3))
# ## mouse.t2 <- stat.t2(mouse.lratio, cl)
# 
# ## Diagnostic plots
# plot.t2(mouse.t2, "Mouse")
# }    
# 
# \keyword{microarray, ttest.}
# 
#*/######################################################################/**


plot.t2 <- function(x, main.title="T plots", low=-5,high=5,...)
{
  par(mfrow=c(2,2),oma=c(1,1,3,1))

  lowt<-x$t < low
  hight<-x$t > high
 
  # 1. t vs. avg. A
  plot(x$A.bar,x$t,xlab="average A",ylab="t",pch=".",
       main="t vs. average A",cex=0.8)
  if(sum.na(lowt)>0)
    points(x$A.bar[lowt],x$t[lowt],pch="*",col=6,cex=1.5)
  if(sum.na(hight)>0)
    points(x$A.bar[hight],x$t[hight],pch="*",col=2,cex=1.5)

   # 2. t_denom vs. avg. A
   plot(x$A.bar,x$Den,xlab="average A",ylab="t denominator",
	pch=".",main="t denominator vs. average A",cex=0.8)
  if(sum.na(lowt)>0)
    points(x$A.bar[lowt],x$Den[lowt],pch="*",col=6,cex=1.5)
  if(sum.na(hight)>0)
    points(x$A.bar[hight],x$Den[hight],pch="*",col=2,cex=1.5)

  # 3. |t_num| vs. avg. A
  plot(x$A.bar,abs(x$Num),xlab="average A",ylab="|t numerator|",
       pch=".",main="|t numerator| vs. average A",cex=0.8)
  if(sum.na(lowt)>0)
    points(x$A.bar[lowt],abs(x$Num[lowt]),pch="*",col=6,cex=1.5)
  if(sum.na(hight)>0)
    points(x$A.bar[hight],abs(x$Num[hight]),pch="*",col=2,cex=1.5)

  # 4. t_denom vs. |t_num|
  plot(abs(x$Num) , x$Den ,xlab="|t numerator|",ylab="t denominator",
       pch=".",main="t denominator vs. |t numerator|",cex=0.8)
  if(sum.na(lowt)>0)
    points(abs(x$Num[lowt]),x$Den[lowt],pch="*",col=6,cex=1.5)
  if(sum.na(hight)>0)
    points(abs(x$Num[hight]),x$Den[hight],pch="*",col=2,cex=1.5)

  par(mfrow=c(1,1))
  mtext(main.title, line=4,cex=1.5)
}

##############################################################
#
# File: Rsingle.R
#
# Created by B. M. Bolstad
# Created on Nov 20,2000 
#
# last modified by bmb
# last modified on Nov 15, 2001
#
# History:
# Nov 15, 2001 Added plot.single.slides function
# Nov 20, 2000 Initial version, created by combining
#              func.churchill.s, func.newton.s
# Jan 18, 2000 Make fixes to complete integration into sma
#
#
###############################################################


#############################################
#
# file :	func.churchill.s
# aim  :	implement methods in
#		Sapir and Churchill
#
# Created:	bmb:	7 June 2000
# Last mod:     bmb:	13 Nov 2000
#
#
# History:	
#	
# 7 June 2000 - initial versions
# 13 Nov 2000 - modifications to allow
#		integration with sma
# 15 Nov 2000 - Continued modification
#
##############################################

##############################################
#
# Internal functions - normalisation functions
# are defunct in sma framework
#
#############################################

	
##############################################
# function to do churchills normalisation.
# input is green and red values
#
# g - green
# r - red
# wts - weights to get robust regression
#
##############################################
	
chur.norm.func <- function(g,r,wts=NA){
# g is green
# r is red

if (is.na(wts)){
    wts <- rep(1,length(g))
}

tmp.fit <- lm(log(r) ~ log(g),weights=wts)
tmp.param <- tmp.fit$coef

#return orthogonal residual
list(param=tmp.param,resid=cos(atan(tmp.param[2]))*resid(tmp.fit))

# misread equation from summary
#tmp.xint <- -tmp.param[1]/tmp.param[2]
#tmp.adj <- log(g) - tmp.xint
#tmp.opp <- fitted(tmp.fit)
#tmp.diff - resid(tmp.fit)*cos(atan(tmp.opp/tmp.adj))

}

##############################################
#
# Wrapper routine to allow input of A and M 
# rather than x,y or g,r to norm function
# basically just converts to x,y then calls
# standard routine
#
# A - inputs
# M - inputs
# weights - to do robust regression
#
###############################################

chur.wrapper.func <- function(A,M,wts=NA){	
	x <- 2^(A - M/2)	
	y <- 2^(A + M/2)
	chur.norm.func(x,y)
}

###############################################
#
# Routine to take input M (log(R/G)) and return
# orthogonal residuals used in Sapir and
# Churchill
#
#
###############################################
	
chur.M.to.e.func <- function(M){
	M/sqrt(2)
}


	
################################################
#
# Function to work if within boundaries basically
# so we can deal with uniform distribution
#
################################################
	
in.func <- function(x,a,b){
(x>=a) & (x<= b)
}


#################################################
#
# Fit Mixture model using EM algorithm
#
# e - orthogonal residuals
# theta - starting parameter estimates
# maxits - maximum iterations for the EM algorithm
# a - lower bounds of uniform distribution
# b - upper bounds of uniform distribution
#
# returns final theta, final posterior probabilities of being different
#	
##################################################
	
chur.em.func <- function(e,theta=c(0.5,1),maxits=50,a=0,b=1)
{
  p <- theta[1]
  s2 <- var(e)
  notdone <- TRUE
  iter <- 1
  while (notdone)
    {
      # E-step
      z <- p*(1/(b-a))*in.func(e,a,b)/(p*in.func(e,a,b)*(1/(b-a))+ (1-p)*(1/sqrt(2*pi*s2)*exp(-e^2/(2*s2))))
      
      # M-step
      p <- sum(z)/length(e)
      s2 <- sum((1-z)*e^2)/(length(e)-sum(z))
#     print(cbind(p,s2))
      iter <- iter +1
      notdone <- (iter < maxits)
    }

  theta[1] <- p
  theta[2] <- s2
  #print(theta)
  list(theta=theta,pp=z)
}

############################################
#
# Churchill mixture model pdf
#
# x - data
# theta - fitted parameters
# a - lower bound for uniform distribution
# b - upper bound for uniform distribution
#
############################################

chur.pdf <- function(x,theta,a=0,b=1)
{
	p <- theta[1]
	s2 <- theta[2]
	(1-p)*(1/sqrt(2*pi*s2)*exp(-x^2/(2*s2))) + p * in.func(x,a,b)*(1/(b-a))
}




#################################################
#
# Function to calculate the values (of M) above 
# which all points have higher posterior probability
# than specified level
#
#
#################################################

givelim <-function(pp,p,s,b,a){
  A <- p*1/(b-a)
  B <- (1-p)*1/sqrt(2*pi*s)
  sqrt(2)*sqrt(-2*log(-A*(pp-1)/(pp*B)))*sqrt(s)
}



###############################################
#
# Only the following functions should be exposed
# to the outside world
#
###############################################
	

#############################################
#
#  Function to perform Churchill Sapir on
#  a set of slides when given 
#
#
#
#############################################

stat.ChurSap <- function(RG,layout,pp=0.95,norm="p", pout=TRUE, image.id=1, ...)
{
  MA <-stat.ma(RG, layout, norm, pout=FALSE)
  stat.ChurSap.ma(MA,pp,pout,image.id,...)
	
}

###############################################
#
# Given we already have M,A (normalised as desired)
# values perform Churchill-Sapir on set of slides
#
################################################
	
stat.ChurSap.ma <- function(MA,pp,pout=TRUE, image.id=1,...)
{
  tmpM <- MA$M[,image.id]
  orthogres <- chur.M.to.e.func(tmpM)
  ind <- !is.na(orthogres)
  pptmp <- rep(NA,length(orthogres))
  orthogres <- orthogres[!is.na(orthogres)]
  A <- min(orthogres)
  B <- max(orthogres)
  #print(orthogres)
  #print(A)
  #print(B)
  em<-chur.em.func(orthogres,a=A,b=B)
#  limits <-givelim(pp,em$theta[1],em$theta[2],B,A)
  if (pout==TRUE){
    plot(MA$A[,image.id],MA$M[,image.id],cex=0.6,xlab="A",ylab="M")
    limits <-givelim(pp,em$theta[1],em$theta[2],B,A)
    abline(limits,0)
    abline(-limits,0)
  } else {
    limits <-givelim(pp,em$theta[1],em$theta[2],B,A)
    pptmp[ind] <- em$pp	
    list(limits=limits,theta=em$theta,pp=pptmp)
  }
}
	


#################################
#
# File:	 func.newton.s
# Aim:	 Implementation of Newtons method
#
# Created ???? June 2000
#
# Last modified:	 Nov 2000
#
# History: 15 Nov 2000 Initial modifications to allow 
#		       integration with sma
#	   18 Nov 2000 Adding deriv information 
#
##################################
	
###########################################################################
# Reading 
# Modify from s.mn2a
###########################################################################

#matt.newton.func <- function(name)
#{
### This based on Newton's normalization.
#  data <- read.table(name,header=T,sep=",") 
#  nspot <- nrow(data)

# Background adjustment (very simple)
#  x <- data[,2]  ## Ch1 green  
#  y <- data[,3]  ## Ch2 red
# Normalization
#  totcy3 <- sum.na( x[x>0] )
#  totcy5 <- sum.na( y[y>0] )
  
#  x <- x/totcy3
#  y <- y/totcy5
  
  
# Rescale to help with underflow problem 10^5 (does not affect shape params)
#  x <- x*100000
#  y <- y*100000
  
#  xx <- x[!is.na(x)]
#  yy <- y[!is.na(y)]

#  list(x=xx, y=yy)
#}

###########################################################################
# newton.func.rotate <- 
# Our modify newton's function.
###########################################################################
###########################################################################
# A) newton.func.rotate <- 
###########################################################################

newton.plot.rotate <- function(A, M, theta,...)
{
## Input A and M
  ind <- is.na(A) | is.na(M) | is.infinite(A) | is.infinite(M)
  A <- A[!ind]
  M <- M[!ind]
  plot(A, M,xlab="A", ylab="M", type="n",...)
  vec1 <- seq(range(A)[1], range(A)[2], length=150)
  vec2 <- seq(range(M)[1], range(M)[2], length=150)
##  theta <- fits[1,]
  logbf <- lod2(A,M,theta=theta)
  bf <- outer(vec1,vec2,"lod2",theta=theta)
  bar <- contour(vec1,vec2,bf,levels=c(0,1,2), save=TRUE, plotit=TRUE, add=TRUE,
		 labex=0, lwd=2 )
  points( A[logbf >=0], M[logbf>=0], cex=.6 , col=2,...)
  points( A[logbf < 0], M[logbf< 0], cex=.6 , col=3,...)
##  box()

}



###########################################################################
# chen.func.rotate <-
###########################################################################
chen.func.rotate <- function(A, M, err=0.01, ...)
{
## reading in A 
## reading in M
  xx <- 2^(A - M/2)
  yy <- 2^(A + M/2)
  tt <- xx/yy
  chat <- sqrt( mean( (tt-1)^2/(1+tt^2) ) )
  tmp01 <- chen.poly(chat, err=err)
 
  abline(h =  log(tmp01[1],2), lty=2, lwd=1.5, ...)
  abline(h =  log(tmp01[2],2), lty=2, lwd=1.5, ...)
}

###########################################################################
# chen.poly (Chen's method)
###########################################################################
chen.poly <- function(cv,err=.01)
        {
        # part of table 2 from Chen et al
        bar <- rbind( c(.979, -2.706, 2.911, -2.805 ),
                      c(.989, 3.082, -2.83, 28.64),
                      c(.9968, -3.496,4.462, -5.002),
                      c( .9648,4.810,-15.161,78.349) )
        if( err==.05 ){
         coef <- bar[1,]
         tmp1 <- cv^3*coef[4] + cv^2*coef[3]+cv*coef[2] + coef[1]
         coef <- bar[2,]
         tmp2 <- cv^3*coef[4] + cv^2*coef[3]+cv*coef[2] + coef[1]
         }
        if( err==.01 )
         {
         coef <- bar[3,]
         tmp1 <- cv^3*coef[4] + cv^2*coef[3]+cv*coef[2] + coef[1]
         coef <- bar[4,]
         tmp2 <- cv^3*coef[4] + cv^2*coef[3]+cv*coef[2] + coef[1]
         }
        return( c(tmp1,tmp2) )
        }


#lod <- function(x,y,theta)
#        {
        # Log_(10) posterior odds
        # x = channel 1 intensity
        # y = channel 2 intensity
        # theta = (aa,a0,nu,pp)
#        aa <- theta[1]; a0 <- theta[2]; x0 <- theta[3]
#         y0 <- x0; z0 <- x0; pp <- theta[4]
#        tmp <- log( pp ) - log(1-pp) +
#                a0*( log(x0) + log(y0) - log(z0) ) +
#                (2*aa+a0)*log(x+y+z0) -
#                (aa+a0)*( log(x+x0) + log(y+y0) ) +
#                2*lgamma(aa+a0) - lgamma(a0) - lgamma(2*aa+a0)
#        return(tmp/2.3)
 #       }

	
###########################################################################
# Calculating lod based on A vs M plot
###########################################################################
lod2 <- function(A, M, theta)
{
# A = log(x*y)/2
# M = log(y/x)
# Log_(10) posterior odds
# x = channel 1 intensity (green)
# y = channel 2 intensity (red)
# theta = (aa,a0,nu,pp)

        x <- 2^(A - M/2)
        y <- 2^(A + M/2)
        aa <- theta[1]
        a0 <- theta[2]
        x0 <- theta[3]
        y0 <- theta[3]
        z0 <- theta[3]
        pp <- theta[4]
        tmp <- log(pp) - log(1 - pp) + a0 * (log(x0) + log(y0) - log(z0)) + (2 * 
                aa + a0) * log(x + y + z0) - (aa + a0) * (log(x + x0) + log(y + 
                y0)) + 2 * lgamma(aa + a0) - lgamma(a0) - lgamma(2 * aa + a0)
        return(tmp/2.3)
}


###########################################################################
# s.em (EM algorithm)
###########################################################################
nploglik <- function(theta,xx=xx,yy=yy,zz=zz)
        {
         # xx,yy are intensities in the two channels; zz=P(b!=c|xx,yy)
         # theta=(aa,a0,nu)
         # (I'll separately optimize pp=P(zz=1); hence npl.. for partial loglik
         aa <- theta[1]; a0<-theta[2]; x0<-theta[3]; y0<- theta[3];
         z0 <- theta[3]
         n <- length(xx)
         # 
         # Complete data loglikelihood
         ll <- (aa-1)*sum(log(xx)+log(yy))  +
                sum(zz)*( 2*(lgamma(aa+a0)-lgamma(aa)-lgamma(a0) ) ) +
                sum(zz)*a0*(log(x0)+log(y0)) +
                (n-sum(zz))*( lgamma(2*aa+a0)-2*lgamma(aa)-lgamma(a0) ) +
                (n-sum(zz))*a0*log(z0) -
                (aa+a0)*sum( zz*( log(x0+xx) + log(y0+yy) ) ) -
                (2*aa+a0)*sum( (1-zz)*( log(z0+xx+yy) ) )
        return(-ll)
        }


func.em <- function(A, M, theta=c(2,1.2,2.7,.4))
{
# EM  algorithm
# input starting values, A and M
# Beta hyperparameter for p
pprior <- 2
# starting value
notdone <- TRUE
iter <- 1
x <- 2^(A - M/2)
y <- 2^(A + M/2)
ind <- is.na(x) | is.na(y) | is.infinite(x) | is.infinite(y)
xx <- x[!ind] 
yy <- y[!ind]
xx <- xx 
yy <- yy

n <- length(xx)
while( notdone )
        {
         aa <- theta[1]; a0<-theta[2]; x0<-theta[3]; y0<- theta[3];
         z0 <- theta[3]; pp <- theta[4]
        # E-step 
        tmp <- log( pp ) - log(1-pp) +
                a0*( log(x0) + log(y0) - log(z0) ) +
                (2*aa+a0)*log(xx+yy+z0) -
                (aa+a0)*( log(xx+x0) + log(yy+y0) ) +
                2*lgamma(aa+a0) - lgamma(a0) - lgamma(2*aa+a0)
        zz <- 1/( 1 + exp(-tmp) )
        # M-step
        #fit <- nlminb( start=c(2,1.2,2.7), objective=nploglik, lower=c(1,0,0),  theta=theta, xx=xx, yy=yy, zz=zz )
	#fit <- optim( par=c(theta[1],theta[2],theta[3]), fn=nploglik,lower=c(1,10^-300,10^-300),gr =nploglikderiv, method="L-BFGS-B",control=list(trace=T),xx=xx, yy=yy, zz=zz )
         fit <- optim( par=c(theta[1],theta[2],theta[3]), fn=nploglik,lower=c(1,10^-300,10^-300),gr =nploglikderiv, method="L-BFGS-B",xx=xx, yy=yy, zz=zz )

         
        # Add a prior on pp
        theta <- c( fit$par,  ( pprior + sum( zz ) )/(2*pprior+n ) )
        #print(round(theta,4) )
        iter <- iter + 1
	 notdone <- (iter < 40)
       } 
theta
}

###########################################################################
# Old function ## without rotation (based on logG vs logR plot
###########################################################################
###########################################################################
# newton.func 
###########################################################################
#
#newton.plot.func <- function(xx, yy, theta,chen=T, chen.err = 0.01)
#{
#  plot( xx, yy, log="xy", pch=".",xlab="Cy3", ylab="Cy5", type="n")
#  vec <- log10(seq(range(c(xx, yy))[1], range(c(xx,yy ))[2], length=150))
##  theta <- fits[1,]
#  logbf <- lod(xx,yy,theta)
#  bf <- outer(10^(vec),10^(vec),"lod",theta=theta)
#  bar <- contour(vec,vec,bf,levels=c(0,1,2), save=T, plotit=T, add=T,
#		 labex=0, lwd=2 )
##  u <- bar$"0"$x
##  v <- bar$"0"$y
##  ind <- is.na(u)
##  u[ind]<- range(vec)
##  v[ind]<- range(vec)
  ## polygon(10^u,10^v,col=4)
#  points( xx[logbf >=0], yy[logbf>=0], cex=.6 , col=2)
#  points( xx[logbf < 0], yy[logbf< 0], cex=.6 , col=3)
##  box()
#  if(chen)
#    chen.func(xx, yy, err=chen.err)
#}


#####################################################################
# chen.func 
#
#
#####################################################################
#chen.func <- function(xx, yy, err=0.01)
#{
#  tt <- xx/yy
#  chat <- sqrt( mean( (tt-1)^2/(1+tt^2) ) )
#  tmp01 <- chen.poly(chat, err=err)
#  # Sandrine changed this: need log10 because plots are in log10 scale
#  abline( log(tmp01[1],10), 1, lty=2, lwd=1.5)
#  abline( log(tmp01[2],10), 1, lty=2, lwd=1.5)
#}


######################################################################
# wrapper function for SMA to perform Newtons method on data
#
#
######################################################################

stat.Newton <- function(RG,layout,norm="p",image.id=1,pout=TRUE){
	MA <-stat.ma(RG, layout, norm, pout=FALSE)
	stat.Newton.ma(MA,image.id,pout)

}


stat.Newton.ma <- function(MA,image.id=1,pout=TRUE){
	M <- MA$M[,image.id]
	A <- MA$A[,image.id]

        ind <- is.na(M) | is.na(A) 
	
	theta <- func.em(A[!ind],M[!ind])

        if (pout == TRUE){
          newton.plot.rotate(A,M,theta)
        } else {
          logodds <- rep(NA,length(M))
          logodds[!ind] <- lod2(A[!ind],M[!ind],theta)
          list(theta=theta,lod=logodds)
        }
      }


#####################################################
#
# nploglikderiv
#
# derivative of loglikelihood function hopefully
# fix broken optim (redundant, found proper scaling)
#
####################################################

	

nploglikderiv <- function(theta,xx=xx,yy=yy,zz=zz){
         # xx,yy are intensities in the two channels; zz=P(b!=c|xx,yy)
         # theta=(aa,a0,nu)
         # (I'll separately optimize pp=P(zz=1); hence npl.. for partial loglik
         aa <- theta[1]; a0<-theta[2]; x0<-theta[3]; y0<- theta[3];
         z0 <- theta[3]
	 nu <- theta[3]
         n <- length(xx)
	daa <- sum(log(xx)+log(yy)) + sum(zz)*(2*digamma(aa+a0)-2*digamma(aa)) +(n-sum(zz))*(2*digamma(2*aa+a0) - 2*digamma(aa)) - sum(zz*(log(x0+xx) + log(y0+yy))) - 2*sum((1-zz)*(log(z0+xx+yy)))
	da0 <- sum(zz)*(2*digamma(aa+a0)-2*digamma(a0)) + sum(zz)*(log(x0)+log(y0)) + (n-sum(zz))*(digamma(2*aa+a0)-digamma(a0)) + (n-sum(zz))*log(z0) - sum(zz*(log(x0+xx)+log(y0+yy))) - sum((1-zz)*log(z0+xx+yy))
	dnu <- 2*sum(zz)*a0/nu + (n-sum(zz))*a0/nu - (aa+a0)*sum(zz*(1/(nu+xx) + 1/(nu+yy))) - (2*aa+a0)*sum((1-zz)/(nu+xx+yy))
	
	c(-daa,-da0,-dnu)
}




chen.plot.rotate <- function(A, M,pout=TRUE){
  ind <- is.na(A) | is.na(M) | is.infinite(A) | is.infinite(M)
  A <- A[!ind]
  M <- M[!ind]
  if (pout==TRUE){
    plot(A, M,cex=0.6,xlab="A", ylab="M")
    chen.func.rotate(A, M, err=0.01, col=4)
    chen.func.rotate(A, M, err=0.05, col=5)
  } else {
      xx <- 2^(A - M/2)
      yy <- 2^(A + M/2)
      tt <- xx/yy
      chat <- sqrt( mean( (tt-1)^2/(1+tt^2) ) )
      tmp01 <- chen.poly(chat, err=0.01)
      tmp05 <- chen.poly(chat, err=0.05)
      list(lower01=log(tmp01[1],2),upper01=log(tmp01[2],2),lower05=log(tmp05[1],2),upper05=log(tmp05[2],2))
  } 
}



stat.Chen <- function(RG,layout,norm="p",image.id=1,pout=TRUE){
  MA <-stat.ma(RG, layout, norm, pout=FALSE)
  stat.Chen.ma(MA,image.id,pout)
}


stat.Chen.ma <- function(MA,image.id,pout=TRUE){
  chen.plot.rotate(MA$A[,image.id],MA$M[,image.id],pout) 
}



plot.single.slide <- function(x,layout,norm="p",image.id=1,...){
  #RG <- x
  MA <- stat.ma(x, layout, norm, pout = FALSE)
  Newton <- stat.Newton.ma(MA,image.id, pout=FALSE)
  ChurSap <- stat.ChurSap.ma(MA,pp=0.95,pout=FALSE,image.id)
  ChurSap2 <- stat.ChurSap.ma(MA,pp=0.99,pout=FALSE,image.id)
  Chen <- stat.Chen.ma(MA,image.id,pout=FALSE)
  newton.plot.rotate(MA$A[,image.id],MA$M[,image.id],Newton$theta,...)
  abline(h=Chen$lower01,col=4,lty=2,lwd=2)
  abline(h=Chen$upper01,col=4,lty=2,lwd=2)
  abline(h=Chen$lower05,col=5,lty=2,lwd=2)
  abline(h=Chen$upper05,col=5,lty=2,lwd=2)
  abline(h=ChurSap$limit,col=6,lty=3,lwd=2)
  abline(h=-ChurSap$limit,col=6,lty=3,lwd=2)
  abline(h=ChurSap2$limit,col=9,lty=3,lwd=2)
  abline(h=-ChurSap2$limit,col=9,lty=3,lwd=2)
}


###########################################################################
# Statistics for Microarray Analysis
# plot.spatial
#
# Date : March 19, 2001
#
# History:
#    March 19, 2001: The spatial plot functions from Rarray.R.
#					
#
# Authors: Sandrine Dudoit and Yee Hwa (Jean) Yang.
##########################################################################

########################################################################/**
# \name{plot.spatial}
#
# \alias{plot.spatial}
# \alias{draw.image.func}
# \alias{spatial.func}
#
# \title{Spatial Representation of Microarray Spot Statistics}
#
#
# \description{Creates an image of shades
# of gray or colors, that represents the values of a statistic for each
# spot on the array. The statistic can be a log intensity ratio, quality
# information such as spot size or shape, or a t-statistic. This function
# can be used to explore whether there are any spatial effects in the data.
# }
#
# \usage{
# plot.spatial(x, layout, crit1=0.15, crit2=crit1, ...)
# }
#
# \arguments{
# \item{x}{a numerical vector. This vector can contain any spot 
# statistics, such
# as log intensity ratios, spot sizes or shapes, or t-statistics.} 
#
# \item{layout}{a list specifying the dimensions of the spot matrix
# and the grid matrix. This can be generated by calling 
# \code{\link{init.grid}}.}
#
# \item{crit1}{the number of values from x to be displayed on the image. 
# If crit1 < 1, the crit1*100\% spots with the largest x values are displayed.
# If crit1 >= 1, the crit1 spots with the largest x values are displayed.}    
#
# \item{crit2}{the number of values from x to be displayed on the
# image. If crit2 < 1, the crit2*100\% spots with the largest x values
# are displayed. If crit2 >= 1, the crit2 spots with the largest x
# values are displayed.}    
# 
# \item{\dots}{graphical parameters may also be supplied as arguments 
# to the function (see \code{\link{par}}).} 
# }
#
# \value{
# An image is created on the current graphics device.
# }
#
# \details{The values that didn't meet the criteria are not shown on the image.
# The image follows the layout of an actual microarray slide with the
# top left corner representing the spot (1,1,1,1).}
#
#
#
# \note{\code{\link{draw.image.func}} and \code{\link{spatial.func}}
#  are called by \code{\link{plot.spatial}} and are not typically 
# used on their own.}
#
# \author{
#  Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
#  Sandrine Dudoit, \email{sandrine@stat.berkeley.edu}}
#
# \seealso{\code{\link{draw.image.func}}, \code{\link{init.grid}}, 
# \code{\link{spatial.func}}, \code{\link{image}}.} 
#
# \examples{
# data(MouseArray)
# # mouse.setup <- init.grid()
# # mouse.data <- init.data() ## see \emph{init.data}
#
# mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# plot.spatial(mouse.lratio$M[,1], mouse.setup) ## default 85% cutoff
#
# # Looking for areas where the spots are not quite circular
# plot.spatial(mouse1[,"shape"], mouse.setup, crit1=0.1)
# }   
#
# \keyword{microarray, spatial.}
#*/#########################################################################

plot.spatial <- function(x, layout, crit1=0.05, crit2=crit1,  ...)
{
  if (crit1 >= 1) crit1 <- crit1 / (length.na(x) - sum(is.na(x)))
  if (crit2 >= 1) crit2 <- crit2 / (length.na(x) - sum(is.na(x)))
  
#  if(crit < 1)
#    tmpind <- x > quantile.na(x, 1-crit)
#  if(crit >=1)
  tmpind <- (x > quantile.na(x, probs=1-crit2)) | (x < quantile.na(x, probs=crit1))

  n <- layout$ngrid.c * layout$ngrid.r * layout$nspot.c * layout$nspot.r

  if (length(as.vector(x)) == n){
    fullm <- x
    fullm[!tmpind] <- NA
  }
  if ((length(as.vector(x)) != n) & (!is.null(names(x)))){
    y <- x[tmpind]; fullm <- rep(NA, n)
    fullm[as.integer(names(y))] <- y
  }
  if ((length(as.vector(x)) != n) & (is.null(names(x)))){
    stop(" Error: Length of vector is different from total number of spots and vector has no row.name.\n")
  }
  draw.image.func(spatial.func(fullm, layout), layout, ...)
}


##########################################################################
# Internal functions called by plot.spatial
##########################################################################

spatial.func <- function(fullm, layout)
{
  gc <- layout$ngrid.c; gr <- layout$ngrid.r
  sc <- layout$nspot.c; sr <- layout$nspot.r
  grid <- split(fullm, rep(1:(gc*gr) , rep(sc*sr, gc*gr)))
  grid1 <- lapply(grid, matrix, nrow=sr, ncol=sc, byrow=TRUE)
  grid2 <- split(unlist(grid1), rep(1:gr, rep(sc*sr*gc, gr)))
  grid3 <- lapply(grid2, matrix, nrow=sr)
  full <- NULL
  for(i in 1:gr){
    full <- rbind(full, grid3[[i]])}
  full
}

draw.image.func <- function (x, layout, axes = FALSE, array.grid = TRUE, label = FALSE, ...) 
{
    gc <- layout$ngrid.c
    gr <- layout$ngrid.r
    sc <- layout$nspot.c
    sr <- layout$nspot.r
    print(summary(as.vector(x)))
    image(1:ncol(x), 1:nrow(x), t(apply(x, 2, rev)), axes = axes, 
        xlab = "", ylab = "", ...)
    if (label) {
        axis(2, ...)
        axis(3, ...)
    }
    box()
    if (array.grid) {
        abline(h = ((gr - 1):1) * (sr) + 0.5)
        abline(v = (1:(gc - 1)) * (sc) + 0.5)
    }
}

############################################################################
#                              End of file
############################################################################
###########################################################################
# Statistics for Microarray Analysis
# T-test
#
# Date : March 19, 2001
#
# History:
#    March 19, 2001: Some of the plot functions from Rarray.R.
#					
#
# Authors: Sandrine Dudoit and Yee Hwa (Jean) Yang.
##########################################################################

##########################################################################
# Test statistics for multiple slides only
##########################################################################

########################################################################/**
# \name{stat.t2}
# 
# \alias{stat.t2}
# \alias{t2stat.func}
# 
# \title{Two-sample t-statistics}
# 
# \description{
# Computes two-sample t-statistics for each gene in a multi-slide
# microarray experiment. 
# }
# 
# \usage{
# stat.t2(X, cl, x.ratio=FALSE, var.equal=TRUE, ...)
# }
# 
# \arguments{
#  \item{X}{if x.ratio=F, X is a list containing two components.  The
#  first component is a matrix of log intensity ratios
#  \eqn{M=\log_2 (R/G)} and the second component is the 
#   average log intensities \eqn{A = log_2 \sqrt{RG}}{A =
#   log_2(R*G)/2}, such as the output 
#   from \code{\link{stat.ma}}. If x.ratio=T, X is a matrix of log
#   expression ratios only.   The rows of X correspond to genes and 
#   columns correspond to different hybridizations, that is different
#   slides.}  
# 
#  \item{cl}{vector of class labels. Must consist of integers 1 and 2.}
# 
#  \item{x.ratio}{logical flag: if TRUE, the matrix X contains only
#  log intensity ratios, if FALSE, X is a list containing two
#  components.  The first component is a matrix of log expression
#  ratios and the second component contains average log
#  intensities A.}
# 
#  \item{var.equal}{logical flag: if TRUE, the variances of the class
#  1 and class 2 parent populations are assumed equal.} 
# 
#  \item{\dots}{other parameters used in \code{\link{t.test}}. }
# }
# 
# \value{
# List containing the following components
# 
#   \item{t}{the two-sample t-statistic for each gene;}
# 
#   \item{Num }{the numerator of the t-statistic for each gene, with
# names attribute "Num";}
# 
#   \item{Denominator}{the denominator of the t-statistic for each gene, with
# names attribute "Den";}
# 
#   \item{n1}{number of class 1 observations used to calculate the
#   t-statistic for each gene;}
# 
#   \item{n2}{number of class 2 observations used to calculate the
#     t-statistics for each gene;}   
# 
#     \item{Average A}{if x.ratio=F, the average across all
#     hybridizations of \eqn{A = log_2 \sqrt{RG}}{A = log_2(R*G)/2}, 
# with names attribute "A.bar", if x.ratio=T, NULL is returned.}
# }
# 
# \references{ D. Freedman, R. Pisani, and
# R. Purves. (1998). Statistics, 3rd ed. NewYork: W.W. Norton.} 
# 
#  
# \note{\code{\link{t2stat.func}} is called by \code{\link{stat.t2}}
# and is not typically used on its own.}         
# 
# \author{
#   Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} \cr
#   Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} 
# }
#   
# \seealso{\code{\link{t2stat.func}}, \code{\link{plot.t2}},
# \code{\link{plot.qq}}, \code{\link{t.test}}.} 
# 
# \examples{
# data(MouseArray)
# ## mouse.setup <- init.grid()
# ## mouse.data <- init.data() ## see \emph{init.data}
# ## mouse.lratio <- stat.ma(mouse.data, mouse.setup)
# cl <- c(rep(1,3), rep(2,3))
# mouse.t2 <- stat.t2(mouse.lratio, cl)
# }
# 
# \keyword{T-test.}
#*/#########################################################################

stat.t2<-function(X, cl, x.ratio=FALSE, var.equal=TRUE,  ...)
{ 
  if(!x.ratio){
    n <- ncol(X)/2
    res<-t(apply(X$M,1,function(z) t2stat.func(z,cl,var.equal, ...)))
    A.bar<-apply(X$A,1,mean.na)
  }
  if(x.ratio){
    n <- ncol(X)
    res<-t(apply(X,1,function(z) t2stat.func(z,cl,var.equal,  ...)))
    A.bar<-NULL
  }
  list(t=res[,"t"],Num=res[,"Num"],Den=res[,"Den"],n1=res[,"n1"],n2=res[,"n2"], A.bar = A.bar)
}


##########################################################################
# Internal Function called by stat.t2
##########################################################################

t2stat.func<-function(x,cl,var.equal=TRUE, ...)
{
  x.ok<-x[!(is.na(x) | is.infinite(x))]
  cl.ok<-cl[!(is.na(x)| is.infinite(x))]
  
  x1<-x.ok[cl.ok==1]
  x2<-x.ok[cl.ok==2]
  n1<-length(x1)
  n2<-length(x2)

  if((n1>2) & (n2>2))
    {
      tmp<-t.test(x1, x2,  var.equal=var.equal, ...)
      tstat<--(tmp$stat)
      num<-tmp$est[2]-tmp$est[1]
      den<-num/tstat
      res<-c(tstat,num,den,n1,n2)
    }
  if((n1<=2) | (n2 <=2))
      res<-c(NA,NA,NA,n1,n2)
  names(res) <- c("t", "Num", "Den", "n1", "n2")
  res
}

############################################################################
#                              End of file
############################################################################
