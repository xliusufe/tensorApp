\name{tmp}
\alias{tmp-function}
\alias{tmp}
\docType{package}
\title{
  Modal Product
}
\description{
  Modal product calculates the product an order \eqn{d} tensor \eqn{S} and a matrix \eqn{M}, that is \eqn{T = S \times _{d0} M} that means \eqn{T_{(d0)} = M \cdot S_{(d0)}}, 
  where \eqn{S_{(d0)}} is the mode-\eqn{d0} unfloding of the tensor \eqn{S}.
}

\usage{tmp(S=NULL, M=NULL, d0=NULL, dims=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{An order \eqn{d} tensor with dimension \code{dims}, where \code{dims} is \eqn{(n_1,\cdots, n_d)}.}
  
  \item{M}{A matrix with \eqn{K} rows and \eqn{n_{d0}} columns. The modal product produces  a new order \eqn{d} tensor with dimension \code{dims} 
           replaced the \code{d0}th dimnesion by \eqn{K}. }
  
  \item{dims}{The dimension of tensor \code{Y}, which is a \eqn{d}-vector \eqn{(n_1,\cdots,n_d)}.}
 
   \item{d0}{\code{d0} is the mode. \eqn{M} multiplies the mode-\code{d0} unfoldings of \code{S}.}  
}

\details{This function gives the modal product of tensor \eqn{S} and matrix \eqn{M}, which is the product of \eqn{M} 
         and the mode-\code{d0} unfolding of tensor \eqn{S}, that is \eqn{S\times _{d0} M=MS_{(d0)}}, 
         where \eqn{S_{(d0)}} is the mode-\code{d0} unfolding of the tensor \eqn{S}. The modal product produce  
         a new order \eqn{d} tensor with dimension \code{dims} replaced the \code{d0} dimnesion by \eqn{K}, where \eqn{K} is the number of rows of matrix \eqn{M}.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Tnew}{Product of \code{S} and \eqn{M}.}
  %\item{...}{ Other options for CompositeQuantile.}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}

\keyword{CP decomposition; HOSVD; Tucker decomposition. }
\seealso{
  TransUnfoldingsT
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
  Y1 <- array(ttu(Y,N,1,dims),dims)
  M <- matrix(1:(4*dims[3]),4)
  
  X1 <- tmp(Y1,M,3,dims)
  X <- ttu(matrix(X1,dims[1]),1,N,dims)
  
  dim(Y1)
  dim(X1)
  print(Y[,1])
  print(X[,1])
  
}