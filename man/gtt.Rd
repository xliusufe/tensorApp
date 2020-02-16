\name{gtt}
\alias{gtt-function}
\alias{gtt}
\docType{package}
\title{Generate a low-rank tensor characterizing the form of Tucker decomposition}
\description{
  This function generates a low-rank tensor characterizing the form of Tucker decomposition with dimension \code{dims} and core tensor \code{S}.
}

\usage{
gtt(dims, S=NULL, d0=NULL, ranks=NULL, seed_id=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

 \item{dims}{The size of the tensor, which is a vector \eqn{(n_1,\cdots,n_d)}.
             \code{dims} must be specified.}
             
  \item{S}{The core tensor. An array with dimension \code{dims}, or a mode-\code{d1}-unfolding of core tensor with size \eqn{n_1\times\cdots\times n_d}.
           Core tensor \code{S} will be given randomly if it is \code{NULL}.}
  
  \item{d0}{\code{d0} is the mode. The output tensor is the mode-\code{d0} unfolding of the tensor. \code{d0} can be NULL (the default) 
            if the output tensor is an array with dimension \code{dims}.} 

  \item{ranks}{The user-specified ranks. It is a vector with length \eqn{d}. If \code{ranks} is \code{NULL} (the default), 
               this function outputs a tensor without low-rank.}
  
  \item{seed_id}{A positive integer, the seed for generating the random numbers. Default is 2.} 
  

}

\details{
  This function generates a low-rank tensor characterizing the form of Tucker decomposition with dimension \code{dims} and core tensor \code{S}.
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
  gtcp, gttsem
}
\examples{
  dims <- c(8,8,10,10,6)
  N <- length(dims)
  ranks <- rep(2,N)
  
  T1 = gtt(dims=dims,d0=1,ranks=ranks)
  T2 <- ttu(T1,1,2,dims)
}
