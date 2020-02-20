# tensorApp
  High-order SVD approximation of a tensor by Tucker or CANDECOMP/PARAFAC (CP) decomposition and rank selection, including sparse PCA by CP decompostion. Transfering a tensor's modal unfoldings to another.
 
  
# Installation

    #install Rtools 3.5 (http://cran.r-project.org/bin/windows/Rtools)
    #install.packages("devtools")
    #install.packages("Rcpp")
    library(devtools)
    install_github("xliusufe/tensorApp")

# Usage

   - [x] [tensorApp-manual](https://github.com/xliusufe/tensorApp/blob/master/inst/tensorApp-manual.pdf) ------------ Details of the usage of the package.
# Example

    library(tensorApp)
    # Example 1 
    # The usage of function "hosvd"
  
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
  
    fit <- hosvd(Y,N,dims,isCP=TRUE)
    Tnew <- fit$Tnew
    ranks1 <- fit$ranks
    lambda <- fit$Tn[[N+1]]
    U1 <- fit$Tn[[1]]
    TNew1 <- ttu(Tnew,N,1,dims)
  
    # Example 2 
    # The usage of function "hosvd_dr"
  
    fit_dr <- hosvd_dr(Y,N,dims,isCP=TRUE)
    Tnew <- fit_dr$Tnew
    ranks1 <- fit_dr$ranks
    lambda <- fit_dr$Tn[[N+1]]
    U1 <- fit_dr$Tn[[1]]
    TNew1 <- ttu(Tnew,N,1,dims)
    
    # Example 3 
    # The usage of function "cpals"
  
    fit_cp <- cpals(Y,N,dims)
    Tnew <- fit_cp$Tnew
    ranks1 <- fit_cp$ranks
    lambda <- fit_cp$Tn[[N+1]]
    U1 <- fit_cp$Tn[[1]]
    TNew2 <- ttu(Tnew,N,1,dims)
	
    # Example 4 
    # The usage of function "ttu"

    T1 <- matrix(1:24,nrow = 4)
    T2 <- ttu(T1,1,2,c(4,3,2))
  
    T0 <- ttu(T2,2,dims=c(4,3,2))  

    # Example 5 
    # The usage of function "cpsym2"  
    dims <- c(8,6,10,6,7)
    N <- length(dims)
    lambda <- seq(6,1,by=-1)
    r1 <- 2
    r2 <- 4
    dr <- 5
    Y = gtcpsem(dims=dims,lambda=lambda,r1=r1,r2=r2,d0=r1,dr=dr)
    fit <- cpsym2(Y,r1=r1,r2=r2,d0=N,dims=dims)
    Tnew <- fit$Tnew
    ranks1 <- fit$ranks
    U1 <- fit$Tn[[r1]]
    U2 <- fit$Tn[[r2]]
    TNew1 <- ttu(Tnew,N,r1,dims)
    TNew2 <- ttu(Tnew,N,r2,dims)  
    
    # Example 6
    library(png)
    bat = readPNG(system.file("data", "bat.png", package="tensorApp"))
    writePNG(bat,target = "bat.png")
    
    Tn = hosvd(bat,dr=20,dims=dim(bat))
    writePNG(Tn$Tnew,target = "batCP.png")
    Tn = hosvd(bat,dr=50,dims=dim(bat),isCP=F,ranks = c(20,20,3))
    writePNG(Tn$Tnew,target = "batTucker.png")
    
    # Example 7
    library(png)
    img = readPNG(system.file("data", "Rlogo.png", package="tensorApp"))
    writePNG(img,target = "Rlogo.png")
    
    Tn = hosvd(img,dr=20,dims=dim(img))
    writePNG(Tn$Tnew,target = "RlogoCP.png")
    Tn = hosvd(img,dr=20,dims=dim(img),isCP=F,ranks = c(20,20,4))
    writePNG(Tn$Tnew,target = "RlogoTucker.png")
    
    # Example 8
    SarsCov2 = readPNG(system.file("data", "SarsCov2.png", package="tensorApp"))
    writePNG(SarsCov2,target = "SarsCov2.png")
    
    Tn = hosvd(SarsCov2,dr=20,dims=dim(SarsCov2))
    writePNG(Tn$Tnew,target = "covid19CP.png")
    Tn = hosvd(SarsCov2,dr=50,dims=dim(SarsCov2),isCP=F,ranks = c(20,20,3))
    writePNG(Tn$Tnew,target = "SarsCov2Tucker.png")  
    
    
    # Example 9
    img = readPNG(system.file("data", "Rlogo.png", package="tensorApp"))
    writePNG(img,target = "Rlogo.png")
    
    Tn = spcacp(img,dr=20,dims=dim(img))
    writePNG(Tn$Tnew,target = "RlogoCP.png")   
  
# Development
The R-package is developed by Xu Liu (liu.xu@sufe.edu.cn).
