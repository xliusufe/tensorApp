
gtcp <- function(dims,lambda=NULL,d0=NULL,dr=10,seed_id=2){

  N = length(dims)
  if(is.null(d0)) d0=0
  if(d0<0||d0>N) stop("d0 must be in the set {1,2,..., N} !")
  set.seed(seed_id)
  if(is.null(lambda))   lambda =  runif(dr,1,3)
  lambda =  sort(lambda,decreasing=TRUE)
  dr = min(c(dr,min(dims),length(lambda)))
  
  D0 = list() 
  for(j in 1:N){
    T1 <- matrix(rnorm(dims[j]^2),nrow = dims[j])
    U <- qr.Q(qr(T1))
    D0[[j]] = U[,1:dr]
  }
  
  dm = prod(dims)
  Dn = matrix(rep(0,dm), nrow = dims[N])
  for(k in 1:dr){
    Uj = D0[[1]][,k]
    for(j in 2:(N-1))   Uj = kronecker(D0[[j]][,k],Uj)
    Dn = Dn + lambda[k]*kronecker(D0[[N]][,k],t(Uj))
  }
  
  if(d0==0) Dn = array(ttu(Dn,N,1,dims),dims)
  else Dn <- ttu(Dn,N,d0,dims)
  
  return (Dn)
}