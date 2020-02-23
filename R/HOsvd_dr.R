
##--------------main by BIC without sparsity----------------------##
hosvd_dr <- function(Y=NULL,d0=NULL,dims=NULL,isCP=TRUE,ranks=NULL,dr=100,D0=NULL,isOrth=FALSE,eps=1e-6,max_step=100,thresh=1e-6){
  
  if(is.null(Y)) stop("Tensor Y must not be NULL !")
  if(is.null(dims)){
    warning("Dimension dims is NULL !") 
    dims = dim(Y)
  }
  N <- length(dims)
  if(is.null(ranks)) ranks = rep(min(dims),N)
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
    D0 = list();
    rmax = max(ranks) 
    for(j in 1:N)    D0[[j]] = diag(dims[j])
    d = rmax^N
    D0[[N+1]] = matrix(runif(d,1,2),rmax)
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
  opts = list(eps=eps,max_step=max_step,N=N,eps1=thresh,max_step1=max_step)
  
  if(isCP){
    if(isOrth) Dnew <- CPTPMorthogon(Y,dd,dr,dims,D0,opts)
    else Dnew <- CPTPM_dr(Y,dd,dr,dims,D0,opts)
    ranks1 = D0[[N+1]]
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
        if(k<dims[j])  ranks1[j] = ranks1[j]+1
        A = D0[[j]]
        D1[[j]] = as.matrix(A[,1:ranks1[j]])
        ds = prod(ranks1)
        A = D0[[N+1]]
        D1[[N+1]] = as.matrix(A[1:ranks1[N],1:(ds/ranks1[N])])
        Dnew <- TuckerALS(Y,dd,dims,ranks1,D1,opts)
        diff = Y - Dnew
        if(sqrt(sum(diff^2)/sum(Y^2))<thresh){
          flag = 1
          break
        }
      }
      if(flag) break
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
