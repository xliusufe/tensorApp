\name{cpsym2}
\alias{cpsym2-function}
\alias{cpsym2}
\docType{package}
\title{
  High-order SVD approximation of a tensor \code{Y} by CP decomposition
}
\description{
  High-order SVD approximation of a semi-symmetric tensor \code{Y} by CANDECOMP/PARAFAC (CP) decomposition with rank preseted or to be selected. 
  The Tensor Power Method is applied. The semi-symmetric tensor means that both mode-\code{r1} and mode-\code{r2} unfoldings are equal, that is, 
  \eqn{Y_{(r1)} = Y_{(r2)}}. 
  For semi-symmetric tensor approximation, the \code{r1}th and \code{r2}th dimensions must be smaller than others.
}

\usage{cpsym2(Y=NULL, r1=1, r2=2, d0=NULL, dims=NULL, dr=10, D0=NULL, isfixr=FALSE, 
                            isOrth=FALSE, eps=1e-4, max_step=50, thresh=1e-6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{An array with dimension \code{dims}, or a \eqn{n_{d0}\times N/n_{d0}} numeric matrix of responses 
           that is the mode \code{d0}-unfolding of tensor in \eqn{\mathcal{R}^{n_1\times\cdots\times n_d}}, where \eqn{N=n_1\times\cdots\times n_d}.}

  \item{r1}{Both \code{r1} and \code{r2} are the user-specified modes, which means 
            that both mode-\code{r1} and mode-\code{r2} unfoldings are equal. Default is \code{r1 = 1} and \code{r2 = 2}.}

  \item{r2}{Both \code{r1} and \code{r2} are the user-specified modes, which means 
            that both mode-\code{r1} and mode-\code{r2} unfoldings are equal. Default is \code{r1 = 1} and \code{r2 = 2}.}
            
  \item{d0}{\code{d0} is the mode. \code{Y} is the mode-\code{d0} unfolding of the tensor. \code{d0} can be NULL (the default) 
            if \code{Y} is an array with dimension \code{dims}.} 
  
  \item{dims}{The size of tensor \code{Y}, which is a \eqn{d}-vector \eqn{(n_1,\cdots,n_d)}. 
             \code{dims} can be NULL (the default) if \code{Y} is an array with dimension \code{dims}. 
             If the length of \code{dims} is 2, it is the ordinary SVD decomposition of a matrix.}
  
  \item{dr}{The user-specified rank for CP decomposition. Default is \code{10}.}

  \item{D0}{A user-specified list of initial matrices of \eqn{U_1,U_2,\cdots,U_d} and core tensor \eqn{S},
            \code{D0=list}\eqn{(U_1=U_1,\cdots,U_d=U_d,S=S)}.
            By default, initial matrices are provided by random.}
            
  \item{isfixr}{A logical value indicating whether the rank is fixed. The rank is  selected automatically if it is \code{FALSE}. Default is \code{FALSE}.}
 
  \item{isOrth}{A logical value indicating whether it outputs orthonognal PCs if CP decomposition is used. Default is \code{FALSE}.}
    
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
    relative change in any coefficient is less than \code{eps}.  Default
    is \code{1e-4}.}
    
  \item{max_step}{Maximum number of iterations.  Default is 50.}
  
  \item{thresh}{Convergence threshhold in the outer loop.  The algorithm iterates until the
             relative change in any coefficient is less than \code{eps}.  Default
             is \code{1e-6}.}    
}

\details{This function gives a \eqn{n_{d0}\times N/n_{d0}} matrix, which is the mode-\code{d0} unfolding, and approximates \code{Y}. 
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Tnew}{Approximation of \code{Y}.}
  
  \item{Tn}{A list of estimated matrices of
              \eqn{U_1,U_2,\cdots,U_d} and core tensor \eqn{S}, \code{Tn=list}\eqn{(U_1=U_1,\cdots,U_d=U_d,S=S)}.}
  
  \item{ranks }{The ranks of estimated tensor \code{Tnew}. It is an integer.}
  %\item{...}{ Other options for CompositeQuantile.}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}

\keyword{CP decomposition; HOSVD; Tucker decomposition. }
\seealso{
  hosvd_dr, tdsym2
}
\examples{ 
  
  dims <- c(6,6,8,7,7)
  N <- length(dims)
  ranks <- rep(2,N)
  dm <- prod(ranks)
  r1 <- 1
  r2 <- 2
  S1 <- matrix(runif(dm,3,7),nrow = ranks[r1])
  S2 <- ttu(S1,r1,r2,ranks)
  S1 <- (S1+S2)/2
  S0 <- ttu(S1,r1,N,ranks)
  
  T1 <- matrix(rnorm(dims[1]*ranks[1]),nrow = dims[1])
  tmp <- qr.Q(qr(T1))
  Uj <- kronecker(tmp,tmp)
  for(k in 3:(N-1)){
    T1 <- matrix(rnorm(dims[k]*ranks[k]),nrow = dims[k])
    tmp <- qr.Q(qr(T1))
    Uj <- kronecker(tmp,Uj)
  }
  T1 <- matrix(rnorm(dims[N]*ranks[N]),nrow = dims[N])
  U <- qr.Q(qr(T1))
  Y <- U\%*\%S0\%*\%t(Uj)
  
  fit <- cpsym2(Y,r1=1,r2=2,d0=N,dims=dims)
  Tnew <- fit$Tnew
  ranks1 <- fit$ranks
  U1 <- fit$Tn[[r1]]
  U2 <- fit$Tn[[r2]]
  TNew1 <- ttu(Tnew,N,r1,dims)
  TNew2 <- ttu(Tnew,N,r2,dims)
}