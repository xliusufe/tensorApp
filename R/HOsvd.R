
assig <- function(n_args){
  cargs <- vector("list", length(n_args))
  for(i in 1:length(n_args)) cargs[[i]] <- 1:n_args[i]
  t(expand.grid(cargs))
}


HOsvd <- function(Y,d0=NULL,dims=NULL,isCP=TRUE,ranks=NULL,dr=20,D0=NULL,eps=1e-6,max_step=100){
  
  if(is.null(Y)) stop("Tensor Y must not be NULL !")
  if(is.null(dims)){
    warning("Dimension dims is NULL !") 
    dims = dim(Y)
  }
  N <- length(dims)
  if(is.null(ranks)) ranks = rep(2,N)
  for(j in 1:N)
    if(ranks[j]>dims[j]) stop(paste("the ",as.character(j),"s rank must be smaller than dimension !",sep = ""))

  if(N<3){
    warning("The input Y is not a tensor with dimension greater than 2! This is an ordinary SVD of a matrix.")
    Ysvd = svd(Y)
    U = Ysvd$u[,1:ranks[1]]
    V = Ysvd$v[,1:ranks[2]]
    d = Ysvd$d
    D = matrix(0,ranks[1],ranks[2])
    diag(D) = d
    Tn = list(U=U,V=V,d=d)
    return (list(Tnew=U%*%D%*%V,
                 Tn=Tn,
                 ranks=ranks
                 )
            )
  } 
  if(is.null(d0)){
    dd = 1
    Y = matrix(Y,dims[1])
  }
  else dd = d0
  if(dd<1||dd>N) stop("d0 must be in the set {1,2,..., N} !")
  
  if(!isCP&is.null(D0)){
    set.seed(1)
    D0 = list() 
    for(j in 1:N){
      U = rbind(diag(ranks[j]), matrix(0,dims[j]-ranks[j],ranks[j]))
      D0[[j]] = U
    }
    S = matrix(runif(prod(ranks),1,2),ranks[N])
    D0[[N+1]] = S
  }
  if(isCP&is.null(D0)){
    set.seed(1)
    D0 = list();
    for(j in 1:N){
      U = rnorm(dims[j])
      D0[[j]] = U/sqrt(sum(U^2))
    }
    D0[[N+1]] = 0
  }
  opts = list(eps=eps,max_step=max_step,N=N,eps1=eps,max_step1=max_step)
  
  if(isCP){
    Dnew <- CPALS(Y,dd,dr,dims,D0,opts)
    ranks <- rep(1,N)
  }
  else Dnew <- TuckerALS(Y,dd,dims,ranks,D0,opts)
  
  if(is.null(d0))  Tnew = array(TransferModalUnfoldingsT(Dnew,N,1,dims),dim=dims)
  else Tnew = TransferModalUnfoldingsT(Dnew,N,d0,dims)
  
  return (list(Tnew=Tnew,
               Tn=D0,
               ranks=ranks
               )
  )
}