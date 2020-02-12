# tensorApp
  High-order SVD approximation by Tucker and CP decomposition and selection of ranks. Transfer a tensor's modal unfoldings to another.
 
  
# Installation

  #install Rtools 3.5 (http://cran.r-project.org/bin/windows/Rtools)
  #install.packages("devtools")
  #install.packages("Rcpp")
  library(devtools)
  install_github("xliusufe/tensorApp")

# Usage

   - [x] [tensorIA-manual](https://github.com/xliusufe/tensorApp/blob/master/inst/tensorApp-manual.pdf) ------------ Details of the usage of the package.
# Example

  library(tensorApp)

  # Example 1 
  # The usage of function "HOsvd()"
  
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
  Y <- U%*%S0%*%t(tmp)
  
  fit <- HOsvd(Y,N,dims,isCP=TRUE)
  Tnew <- fit$Tnew
  ranks1 <- fit$ranks
  TNew1 <- TransUnfoldingsT(Tnew,N,1,dims)
  
  # Example 2 
  # The usage of function "HOsvd_dr()"
  
  fit_dr <- HOsvd_dr(Y,N,dims,isCP=TRUE)
  Tnew <- fit_dr$Tnew
  ranks1 <- fit_dr$ranks
  TNew1 <- TransUnfoldingsT(Tnew,N,1,dims)

  # Example 3 
  # The usage of function "TransUnfoldingsT()"

  T1 <- matrix(1:24,nrow = 4)
  T2 <- TransUnfoldingsT(T1,1,2,c(4,3,2))
  
  T0 <- TransUnfoldingsT(T2,2,dims=c(4,3,2))  
  
# Development
The R-package is developed by Xu Liu (liu.xu@sufe.edu.cn).
