\name{HOsvd}
\alias{HOsvd-function}
\alias{HOsvd}
\docType{package}
\title{
  High-order SVD approximation by Tucker and CP decomposition
}
\description{
  High-order SVD approximation by Tucker decomposition or CANDECOMP/PARAFAC (CP) decomposition with preset rank.
}

\usage{HOsvd(Y, d0=NULL, dims=NULL, isCP=TRUE, ranks=NULL, dr=20, 
                                                D0=NULL, eps=1e-6, max_step=100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{An array with dimension \code{dims}, or a \eqn{n_{d0}\times N/n_{d0}} numeric matrix of responses 
           that is the mode \code{d0}-unfolding of tensor in \eqn{\mathcal{R}^{n_1\times\cdots\times n_d}}, where \eqn{N=n_1\times\cdots\times n_d}.}

  \item{d0}{\code{d0} is the mode which unfoldings \code{Y} is. \code{d0} can be NULL (the default) if \code{Y} is an array with dimension \code{dims}.} 
  
  \item{dims}{The size of tensor \code{Y}, which is a \eqn{d}-vector \eqn{(n_1,\cdots,n_d)}. 
             \code{dims} can be NULL (the default) if \code{Y} is an array with dimension \code{dims}. 
             If the length of \code{dims} is 2, it is the ordinary SVD decomposition of a matrix.}
  
  \item{isCP}{A logical value indicating whether CP decomposition will be used. Default is \code{TRUE}.}
  
  \item{ranks}{The user-specified ranks. It is a vector with length \eqn{d}. Default is \eqn{(2,\cdots,2)}.}
  
  \item{dr}{The user-specified rank for CP decomposition. It is useless if Tucker decomposition is used. Default is 20.}
  
  \item{D0}{A user-specified list of initial matrices of \eqn{U_1,U_2,\cdots,U_d} and core tensor \eqn{S},
            \code{D0=list}\eqn{(U_1=U_1,\cdots,U_d=U_d,S=S)}.
            By default, initial matrices are provided by random.}
    
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
    relative change in any coefficient is less than \code{eps}.  Default
    is \code{1e-6}.}
    
  \item{max_step}{Maximum number of iterations.  Default is 100.}
}

\details{This function gives a \eqn{n_{d0}\times N/n_{d0}} matrix, which is the mode-\code{d0} unfolding, and approximates \code{Y}. 
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Tnew}{Approximation of \code{Y}.}
  
  \item{Tn}{A list of estimated matrices of
              \eqn{U_1,U_2,\cdots,U_d} and core tensor \eqn{S}, \code{Tn=list}\eqn{(U_1=U_1,\cdots,U_d=U_d,S=S)}.}
  
  \item{ranks }{The ranks of estimated tensor \code{Tnew}.}
  %\item{...}{ Other options for CompositeQuantile.}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}

\keyword{CP decomposition; HOSVD; Tucker decomposition. }
\seealso{
  HOsvd_dr
}
\examples{ 
  
  dims <- c(8,8,10,10,6)
  N <- length(dims)
  ranks <- rep(2,N)
  S0 = matrix(runif(prod(ranks),3,7),ranks[N])
  T1 <- matrix(rnorm(dims[1]*ranks[1]),nrow = dims[1])
  tmp <- qr.Q(qr(T1))
  for(k in 2:(N-1)){
    T1 <- matrix(rnorm(dims[k]*ranks[k]),nrow = dims[k])
    tmp <- kronecker(qr.Q(qr(T1)),tmp)
  }
  T1 <- matrix(rnorm(dims[N]*ranks[N]),nrow = dims[N])
  U = qr.Q(qr(T1))
  Y <- U\%*\%S0\%*\%t(tmp)
  
  fit <- HOsvd(Y,N,dims,isCP=TRUE)
  Tnew <- fit$Tnew
  ranks1 <- fit$ranks
  TNew1 <- TransUnfoldingsT(Tnew,N,1,dims)
}