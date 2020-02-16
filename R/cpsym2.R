
cpsym2 <- function(Y=NULL,r1=1,r2=2,d0=NULL,dims=NULL,dr=10,D0=NULL,isfixr=TRUE,isOrth=FALSE,eps=1e-4,max_step=50,thresh=1e-6){
  
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
    for(j in 1:N){
      if(j!=r1&&j!=r2){
        U = runif(dims[j],1,3)
        D0[[j]] = U/sqrt(sum(U^2))
      }
    }
    U = runif(dims[j],1,3)
    D0[[r1]] = U/sqrt(sum(U^2))
    D0[[r2]] = D0[[r1]]
    S = 0
    D0[[N+1]] = S
  }
  if(isOrth)   Dnew <- CPTPMsym2Orth(Y,dr,r1-1,r2-1,dims,D0,opts)
  else  Dnew <- CPTPMsym2(Y,dr,r1-1,r2-1,dims,D0,opts)
  
  ranks <- D0[[N+1]]
  
  if(is.null(d0))  Tnew = array(TransferModalUnfoldingsT(Dnew,N,1,dims),dim=dims)
  else Tnew = TransferModalUnfoldingsT(Dnew,N,d0,dims)
  
  return (list(Tnew=Tnew,
               Tn=D0,
               ranks=ranks
               )
  )
}