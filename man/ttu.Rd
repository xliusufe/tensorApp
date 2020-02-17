\name{ttu}
\alias{ttu-function}
\alias{ttu}
\docType{package}
\title{Transfer a tensor's modal unfolding to another.}
\description{
  Transfer a tensor's modal unfolding to another. 
}

\usage{
ttu(S=NULL, d1=NULL, d2=0, dims=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{An array with dimension \code{dims}, or a mode-\code{d1}-unfolding of a tensor with size \eqn{n_1\times\cdots\times n_d}.}
  
  \item{d1}{An integer, the mode of unfolding \eqn{S_{(d_1)}}.  \code{d1} can be NULL (the default) if \code{S} is an array with dimension \code{dims}.}
  
  \item{d2}{An integer, the mode of output unfolding \eqn{S_{(d_2)}}. It transfers \code{S} to an array with dimension \code{dims} if \code{d2=0}.
            The default is \code{0}.}
  
 \item{dims}{The size of tensor \eqn{S}, which is a vector \eqn{(n_1,\cdots,n_d)}.
             \code{dims} can be NULL (the default) if \code{S} is an array with dimension \code{dims}.}
}

\details{
  This function transfers an input mode-\code{d1}-unfolding \eqn{S_{(d_1)}} 
  to mode-\code{d2}-unfolding \eqn{S_{(d_2)}}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Td2}{the output mode-\code{d2}-unfolding, \eqn{S_{(d_2)}}.}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}

\keyword{ CP decomposition; HOSVD; Tucker decomposition}
\examples{
  T1 <- matrix(1:24,nrow = 4) # A tensor unfolding with size 4*6
  T2 <- ttu(T1,1,2,c(4,3,2))
  
  T0 <- ttu(T2,2,dims=c(4,3,2))
}
