\name{hosvd}
\alias{hosvd-function}
\alias{hosvd}
\docType{package}
\title{
  High-order SVD approximation of a tensor \code{Y} by Tucker or CP decomposition
}
\description{
  High-order SVD approximation of a tensor \code{Y} by Tucker decomposition or CANDECOMP/PARAFAC (CP) decomposition with preseted rank. 
  Alternating Least Squares algorithm is applied to Tucker decomposition, and Tensor Power Method is applied to CP decomposition. 
}

\usage{hosvd(Y=NULL, d0=NULL, dims=NULL, isCP=TRUE, ranks=NULL, dr=20, 
                                  D0=NULL,isOrth=FALSE, eps=1e-6, max_step=100)
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
  
  \item{isCP}{A logical value indicating whether CP decomposition will be used. Default is \code{TRUE}.}
  
  \item{ranks}{The user-specified ranks. It is a vector with length \eqn{d}. Default is \eqn{(2,\cdots,2)}.}
  
  \item{dr}{The user-specified rank for CP decomposition. It is useless if Tucker decomposition is used. Default is 20.}
  
  \item{D0}{A user-specified list of initial matrices of \eqn{U_1,U_2,\cdots,U_d} and the mode-\eqn{d0} unfolding \eqn{S_{(d0)}} of the core tensor \eqn{S},
            \code{D0=list}\eqn{(U_1=U_1,\cdots,U_d=U_d,S=S_{(d0)})}.
            By default, initial matrices are provided randomly.}
    
  \item{isOrth}{A logical value indicating whether it outputs orthonognal PCs if CP decomposition is used. Default is \code{FALSE}.}
  
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
  
  \item{Tn}{A list of estimated matrices of \eqn{U_1,U_2,\cdots,U_d} and the mode-\eqn{d0} unfolding \eqn{S_{(d0)}} of the core tensor \eqn{S}, 
            \code{Tn=list}\eqn{(U_1=U_1,\cdots,U_d=U_d,S=S_{(d0)})}.}
  
  \item{ranks }{The ranks of estimated tensor \code{Tnew}. It is a vector with the same length as \code{dims} if Tucker decomposition is used, 
                or an integer if CP decomposition is used.}
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
  
  # Example 1
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
  
  fit <- hosvd(Y,N,dims,isCP=TRUE)
  Tnew <- fit$Tnew
  ranks1 <- fit$ranks
  lambda <- fit$Tn[[N+1]]
  U1 <- fit$Tn[[1]]
  TNew1 <- ttu(Tnew,N,1,dims)
  
  
  # Example 2
  library(png)
  bat = readPNG(system.file("data", "bat.png", package="tensorApp"))
  writePNG(bat,target = "bat.png")
  
  Tn = hosvd(bat,dr=20,dims=dim(bat))
  writePNG(Tn$Tnew,target = "batCP.png")
  Tn = hosvd(bat,dr=50,dims=dim(bat),isCP=F,ranks = c(20,20,3))
  writePNG(Tn$Tnew,target = "batTucker.png")
  
  # Example 3
  img = readPNG(system.file("data", "Rlogo.png", package="tensorApp"))
  writePNG(img,target = "Rlogo.png")
  
  Tn = hosvd(img,dr=20,dims=dim(img))
  writePNG(Tn$Tnew,target = "RlogoCP.png")
  Tn = hosvd(img,dr=20,dims=dim(img),isCP=F,ranks = c(20,20,4))
  writePNG(Tn$Tnew,target = "RlogoTucker.png")
  
  
  # Example 4
  SarsCov2 = readPNG(system.file("data", "SarsCov2.png", package="tensorApp"))
  writePNG(SarsCov2,target = "SarsCov2.png")

  Tn = hosvd(SarsCov2,dr=20,dims=dim(SarsCov2))
  writePNG(Tn$Tnew,target = "SarsCov2CP.png")
  Tn = hosvd(SarsCov2,dr=50,dims=dim(SarsCov2),isCP=F,ranks = c(20,20,3))
  writePNG(Tn$Tnew,target = "SarsCov2Tucker.png")
}