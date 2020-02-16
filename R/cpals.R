
cpals <- function(Y=NULL,d0=NULL,dims=NULL,dr=NULL,D0=NULL,eps=1e-4,max_step=50,thresh=1e-6){
  
  if(is.null(Y)) stop("Tensor Y must not be NULL !")
  if(is.null(dims)){
    warning("Dimension dims is NULL !") 
    dims = dim(Y)
  }
  N <- length(dims)
  if(N<3){
    warning("The input Y is not a tensor with dimension greater than 2! This is an ordinary SVD of a matrix.")
    Ysvd = svd(Y)
    U = Ysvd$u
    V = Ysvd$v
    d = Ysvd$d
    D = matrix(0,nrow(Y),ncol(Y))
    diag(D) = d
    Tn = list(U=U,V=V,d=d)
    return (list(Tnew=U%*%D%*%V,
                 Tn=Tn,
                 ranks=dim(Y)
                 )
            )
  } 
  if(is.null(d0)){
    dd = 1
    Y = matrix(Y,dims[1])
  }
  else dd = d0
  if(dd<1||dd>N) stop("d0 must be in the set {1,2,..., N} !")
  
  opts = list(eps=eps,max_step=max_step,N=N)
  
  if(!is.null(dr)){
    if(is.null(D0)){
      set.seed(1)
      D0 = list() 
      for(j in 1:N){
        U = matrix(runif(dims[j]*dr,1,3),dims[j],dr)
        U[1,] = 1
        D0[[j]] = U
      }
      S = 0
      D0[[N+1]] = S
    }
    Dnew <- CPALS(Y,dd,dr,dims,D0,opts)
  }
  else{
    dm = 100
    if(is.null(D0)){
      set.seed(1)
      D0 = list() 
      for(j in 1:N){
        U = matrix(runif(dims[j]*dm,1,3),dims[j],dm)
        U[1,] = 1
        D0[[j]] = U
      }
      S = 0
      D0[[N+1]] = S
    }
    D1 = list()
    for(k in 2:dm){
      for(j in 1:N){
        U = D0[[j]]
        D1[[j]] = U[,1:k]
      }
      Dnew <- CPALS(Y,dd,k,dims,D1,opts)
      diff = Y - Dnew
      if(sqrt(sum(diff^2)/sum(Y^2))<thresh) break
    }
  }
  
  if(is.null(d0))  Tnew = array(TransferModalUnfoldingsT(Dnew,N,1,dims),dim=dims)
  else Tnew = TransferModalUnfoldingsT(Dnew,N,d0,dims)
  
  return (list(Tnew=Tnew,
               Tn=D1,
               ranks=k
               )
  )
}