\name{spcacp}
\alias{spcacp-function}
\alias{spcacp}
\docType{package}
\title{
  High-order SVD approximation of a tensor \code{Y} by sparse CP decomposition
}
\description{
  High-order SVD approximation of a tensor \code{Y} by sparse CANDECOMP/PARAFAC (CP) decomposition with rank preseted or to be selected. 
  The Alternating Least Squares (als) algorithm is applied.
}

\usage{spcacp(Y=NULL, d0=NULL, dims=NULL, nactive=NULL, dr=10, criteria="BIC", 
              penalty="LASSO", D0=NULL, lambda=NULL, nlam=50, lam_min=1e-4,
              gamma=2, alpha=1, eps=1e-4, max_step=20)
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
  
  \item{nactive}{It is an integer vector, which are dimensional indeces without requirement of penalization. 
                Default is NULL, which means all dimensions will be penalized.}
  
  \item{dr}{The user-specified rank for CP decomposition. Default is \code{10}.}

  \item{D0}{A user-specified list of initial matrices of \eqn{U_1,U_2,\cdots,U_d} and core tensor \eqn{S},
            \code{D0=list}\eqn{(U_1=U_1,\cdots,U_d=U_d,S=S)}.
            By default, initial matrices are provided by random.}
    
    
  \item{criteria}{The criteria to be applied to select tuning parameters.  Either \code{BIC} (the default), \code{AIC}, 
                or \code{GCV}.}

  \item{penalty}{The penalty to be applied to the model. Either "LASSO" (the default),  
                 "SCAD", or "MCP".}
        
  \item{lambda}{A user-specified sequence of lambda values.  By default,
        a sequence of values of length \code{nlam} is computed, equally
        spaced on the log scale.}

  \item{nlam}{The number of lambda values. Default is 50.}
  
  \item{lam_min}{The smallest value for lambda, as a fraction of
                 lambda.max.  Default is 1e-4.}

  \item{gamma}{The tuning parameter of the MCP/SCAD penalty (see details).}
  
  \item{alpha}{Tuning parameter for the Mnet estimator which controls
               the relative contributions from the LASSO, MCP/SCAD penalty and the ridge,
               or L2 penalty.  \code{alpha=1} is equivalent to LASSO, MCP/SCAD penalty,
               while \code{alpha=0} would be equivalent to ridge regression.
               However, \code{alpha=0} is not supported; \code{alpha} may be
               arbitrarily small, but not exactly 0.}
               
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
    relative change in any coefficient is less than \code{eps}.  Default
    is \code{1e-4}.}
    
  \item{max_step}{Maximum number of iterations.  Default is 20.}
  
}

\details{This function gives a \eqn{n_{d0}\times N/n_{d0}} matrix, which is the mode-\code{d0} unfolding, and approximates sparse tensor \code{Y}. 
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Tnew}{Approximation of \code{Y}.}
  
  \item{Tn}{A list of estimated matrices of
              \eqn{U_1,U_2,\cdots,U_d} and core tensor \eqn{S}, \code{Tn=list}\eqn{(U_1=U_1,\cdots,U_d=U_d,S=S)}.}
  
  \item{ranks }{The ranks of estimated tensor \code{Tnew}. It is an integer.}
  
  \item{rss}{Residual sum of squares (RSS).}
  
  \item{lambda}{The sequence of regularization parameter values in the path.}
  
  \item{selectedID}{The index of \code{lambda} corresponding to
                    \code{lambda_opt}.}
  
  \item{lambda_opt}{The value of \code{lambda} with the minimum
                    \code{BIC} value.}
  
  \item{df}{Degrees of freedom.}  
  %\item{...}{ Other options for CompositeQuantile.}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}

\keyword{CP decomposition; HOSVD; Sparse PCA; Tucker decomposition. }
\seealso{
  pcals
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
  
  fit <- spcacp(Y,N,dims)
  Tnew <- fit$Tnew
  ranks1 <- fit$ranks
  U1 <- fit$Tn[[1]]
  TNew1 <- ttu(Tnew,N,1,dims)
  
  
  # Example 2
  img = readPNG(system.file("img", "Rlogo.png", package="tensorApp"))
  writePNG(img,target = "Rlogo.png")
  
  Tn = spcacp(img,dr=20,dims=dim(img))
  writePNG(Tn$Tnew,target = "RlogoCP.png")
  
}