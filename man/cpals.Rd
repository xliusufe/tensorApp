\name{cpals}
\alias{cpals-function}
\alias{cpals}
\docType{package}
\title{
  High-order SVD approximation of a tensor \code{Y} by CP decomposition
}
\description{
  High-order SVD approximation of a tensor \code{Y} by CANDECOMP/PARAFAC (CP) decomposition with preset rank or rank to be selected. 
  The Alternating Least Squares (als) algorithm is applied.
}

\usage{cpals(Y=NULL, d0=NULL, dims=NULL, dr=10, isfixr=FALSE, 
                              D0=NULL, eps=1e-4, max_step=50, thresh=1e-6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{An array with dimension \code{dims}, or a \eqn{n_{d0}\times N/n_{d0}} numeric matrix of responses 
           that is the mode \code{d0}-unfolding of tensor in \eqn{\mathcal{R}^{n_1\times\cdots\times n_d}}, where \eqn{N=n_1\times\cdots\times n_d}.}

  \item{d0}{\code{d0} is the mode. \code{Y} is the mode-\code{d0} unfolding of the tensor. \code{d0} can be NULL (the default) 
            if \code{Y} is an array with dimension \code{dims}.} 
  
  \item{dims}{The size of tensor \code{Y}, which is a \eqn{d}-vector \eqn{(n_1,\cdots,n_d)}. 
             \code{dims} can be NULL (the default) if \code{Y} is an array with dimension \code{dims}. 
             If the length of \code{dims} is 2, it is the ordinary SVD decomposition of a matrix.}
  
  \item{dr}{The user-specified rank for CP decomposition. Default is \code{10}.}

 \item{isfixr}{A logical value indicating whether the rank is fixed. The rank is  selected automatically if it is \code{FALSE}. Default is \code{FALSE}.}
 
  \item{D0}{A user-specified list of initial matrices of \eqn{U_1,U_2,\cdots,U_d} and core tensor \eqn{S},
            \code{D0=list}\eqn{(U_1=U_1,\cdots,U_d=U_d,S=S)}.
            By default, initial matrices are provided by random.}
    
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
  hosvd_dr
}
\examples{ 
  
  dims <- c(8,8,10,10,6)
  N <- length(dims)
  ranks <- rep(2,N)
  S0 <- matrix(runif(prod(ranks),3,7),ranks[N])
  T1 <- matrix(rnorm(dims[1]*ranks[1]),nrow = dims[1])
  tmp <- qr.Q(qr(T1))
  for(k in 2:(N-1)){
    T1 <- matrix(rnorm(dims[k]*ranks[k]),nrow = dims[k])
    tmp <- kronecker(qr.Q(qr(T1)),tmp)
  }
  T1 <- matrix(rnorm(dims[N]*ranks[N]),nrow = dims[N])
  U <- qr.Q(qr(T1))
  Y <- U\%*\%S0\%*\%t(tmp)
  
  fit <- cpals(Y,N,dims)
  Tnew <- fit$Tnew
  ranks1 <- fit$ranks
  U1 <- fit$Tn[[1]]
  TNew1 <- ttu(Tnew,N,1,dims)
  
}