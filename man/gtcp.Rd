\name{gtcp}
\alias{gtcp-function}
\alias{gtcp}
\docType{package}
\title{Generate a low-rank tensor characterizing the form of CP decomposition}
\description{
  This function generates a low-rank tensor characterizing the form of CP decomposition with dimension \code{dims} and factors \code{lambda}.
}

\usage{
gtcp(dims, lambda=NULL, d0=NULL, dr=NULL, seed_id=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

 \item{dims}{The size of the tensor, which is a vector \eqn{(n_1,\cdots,n_d)}.
             \code{dims} must be specified.}
             
  \item{lambda}{The factors of CP decomposition. It is an vector. Factors \code{lambda} will be given randomly if it is \code{NULL}.}
  
  \item{d0}{\code{d0} is the mode. The output tensor is the mode-\code{d0} unfolding of the tensor. \code{d0} can be NULL (the default) 
            if the output tensor is an array with dimension \code{dims}.} 

  \item{dr}{The user-specified rank. Default is \code{10}.}
  
  \item{seed_id}{A positive integer, the seed for generating the random numbers. Default is 2.} 
  

}

\details{
  This function generates a low-rank tensor characterizing the form of CP decomposition with dimension \code{dims} and factors \code{lambda}.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Dn}{the output mode-\code{d0}-unfolding, \eqn{D_{(d_0)}}. Or an array with dimesion \code{dims} if \code{d0} is \code{NULL}.}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}

\keyword{ CP decomposition; HOSVD; Tucker decomposition}
\seealso{
  gtcpsem, gtt
}
\examples{
  dims <- c(8,8,10,10,6)
  N <- length(dims)
  lambda <- seq(6,1,by=-1)
  dr <- 5
  
  T1 <- gtcp(dims=dims,lambda=lambda,d0=1,dr=dr)
  T2 <- ttu(T1,1,2,dims)
}
