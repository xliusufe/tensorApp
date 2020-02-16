
gtt <- function(dims,S=NULL,d0=NULL,ranks=NULL,seed_id=2){

  N = length(dims)
  if(is.null(d0)) d0=0
  if(d0<0||d0>N) stop("d0 must be in the set {1,2,..., N} !")
  set.seed(seed_id)
  if(is.null(ranks)){
    dm = prod(dims)
    S1 = matrix(runif(dm,3,7),dims[1])
    if(d0==0)     Dn = array(S1,dims)
    else Dn = ttu(S1,1,d0,dims)
    return(Dn)
  }
  if(N!=length(ranks)) stop("Both lengths of dimensions and ranks must be equal !")
  
  dm = prod(ranks)
  if(is.null(S))  S0 = matrix(runif(dm,3,7),ranks[N])
  else{
    rm = prod(dim(S))
    if(rm!=dm) stop("Core tensor S dismatch the ranks!")
    S0 = matrix(S,ranks[N])
  }
  
  T1 <- matrix(rnorm(dims[1]*ranks[1]),nrow = dims[1])
  Uj <- qr.Q(qr(T1))
  for(k in 2:(N-1)){
    T1 <- matrix(rnorm(dims[k]*ranks[k]),nrow = dims[k])
    tmp <- qr.Q(qr(T1))
    Uj <- kronecker(tmp,Uj)
  }
  T1 <- matrix(rnorm(dims[N]*ranks[N]),nrow = dims[N])
  Un = qr.Q(qr(T1))
  Dn <- Un%*%S0%*%t(Uj)
  if(d0==0) Dn = array(ttu(Dn,N,1,dims),dims)
  else Dn <- ttu(Dn,N,d0,dims)
  
  return (Dn)
}