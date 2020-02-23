
tdsym2 <- function(Y=NULL,r1=1,r2=2,d0=NULL,dims=NULL,ranks=NULL,dr=10,D0=NULL,isfixr=TRUE,eps=1e-4,max_step=50,thresh=1e-6){

  if(r1==r2){
    warning("Two modes are the same ! The tensor is approximated without symmetry.")
    fit = hosvd_dr(Y,d0,dims,ranks=ranks,dr=dr,D0=D0,isOrth=isOrth)
    return(fit)
  }
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
  
  opts = list(eps=eps,max_step=max_step,N=N,eps1=thresh,isfixr=as.numeric(isfixr))
  
  if(dd!=N) Y = TransferModalUnfoldingsT(Y,dd,N,dims)
  if(!isfixr)    dr = max(dr,100)
  if(is.null(D0)){
    set.seed(1)
    D0 = list() 
    rmax = max(ranks) 
    for(j in 1:N)
      if(j!=r1&&j!=r2)    D0[[j]] = diag(dims[j])
    U = diag(dims[r1])
    D0[[r1]] = U
    D0[[r2]] = U
    
    dm = prod(ranks)
    S1 = matrix(runif(dm,3,7),ranks[r1])
    S0 = ttu(gtsem0(S1,r1,r2,ranks),r1,N,ranks)
    D0[[N+1]] = S0
  }
  if(isfixr){
    Dnew <- TuckerALSsym2(Y,r1-1,r2-1,dims,ranks,D0,opts)
    ranks1 = ranks
    D1 = D0
  }
  else{
    ranks1 = rep(1,N)
    D1 = list()
    for(j in 1:N){
      A = D0[[j]]
      D1[[j]] = as.matrix(A[,1])
    }
    D1[[N+1]] = 1
    flag = 0
    for(k in 1:dr){
      for(j in 1:N){
        if(j!=r1&j!=r2){
          if(k<dims[j])  ranks1[j] = ranks1[j]+1
          A = D0[[j]]
          D1[[j]] = as.matrix(A[,1:ranks1[j]])
          ds = prod(ranks1)
          A = D0[[N+1]]
          D1[[N+1]] = as.matrix(A[1:ranks1[N],1:(ds/ranks1[N])])
          Dnew <- TuckerALSsym2(Y,r1-1,r2-1,dims,ranks1,D1,opts)
          diff = Y - Dnew
          if(sqrt(sum(diff^2)/(sum(Y^2)+1))<thresh){
            flag = 1
            break
          }
        }
      }
      if(flag) break
      
      j=r1
      if(k<dims[j]){
        ranks1[j] = ranks1[j]+1
        ranks1[r2] = ranks1[r2]+1
      }
      A = D0[[j]]
      D1[[j]] = as.matrix(A[,1:ranks1[j]])
      ds = prod(ranks1)
      A = D0[[N+1]]
      D1[[N+1]] = as.matrix(A[1:ranks1[N],1:(ds/ranks1[N])])
      Dnew <- TuckerALSsym2(Y,r1-1,r2-1,dims,ranks1,D1,opts)
      diff = Y - Dnew
      if(sqrt(sum(diff^2)/(sum(Y^2)+1))<thresh)   break
    }
  }
  
  if(is.null(d0))  Tnew = array(TransferModalUnfoldingsT(Dnew,N,1,dims),dim=dims)
  else Tnew = TransferModalUnfoldingsT(Dnew,N,d0,dims)
  
  return (list(Tnew=Tnew,
               Tn=D1,
               ranks=ranks1
               )
  )
}