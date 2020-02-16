
gttsem <- function(dims,S=NULL,r1=1,r2=2,d0=NULL,ranks=NULL,seed_id=2){
  
  if(dims[r1]!=dims[r2]) stop("Both the r1 and r2 dimensions must be equal !")
  N = length(dims)
  if(is.null(d0)) d0=0
  if(d0<0||d0>N) stop("d0 must be in the set {1,2,..., N} !")
  if(r1>r2){
    r = r1
    r1 = r2
    r2 = r    
  }
  set.seed(seed_id)
  if(is.null(ranks)){
    dm = prod(dims)
    S1 = matrix(runif(dm,3,7),dims[r1])
    S1 = gtsem0(S1,r1,r2,dims)
    if(d0==0)   Dn = array(ttu(S1,r1,1,dims),dims)
    else Dn = ttu(S1,r1,d0,dims)
    return(Dn)
  }
  if(N!=length(ranks)) stop("Both lengths of dimensions and ranks must be equal !")
  if(ranks[r1]!=ranks[r2]) stop("Both the r1 and r2 ranks must be equal !")
  if(r2-r1>2) stop("The difference of r1 and r2 must be smaller than 3 !")
  
  dm = prod(ranks)
  if(is.null(S))  S1 = matrix(runif(dm,3,7),ranks[r1])
  else{
    rm = prod(dim(S))
    if(rm!=dm) stop("Core tensor S dismatch the ranks!")
    S1 = matrix(S,ranks[r1])
  }
  S0 = ttu(gtsem0(S1,r1,r2,ranks),r1,N,ranks)

  T1 <- matrix(rnorm(dims[1]*ranks[1]),nrow = dims[1])
  Uj <- qr.Q(qr(T1))
  if(r1==1) Ur1 = Uj
  for(k in 2:(N-1)){
    if(k==r2) tmp = Ur1
    else{
      T1 <- matrix(rnorm(dims[k]*ranks[k]),nrow = dims[k])
      tmp <- qr.Q(qr(T1))
      if(k==r1) Ur1 = tmp
    }
    Uj <- kronecker(tmp,Uj)
  }
  if(N==r2) Un = Ur1
  else{
    T1 <- matrix(rnorm(dims[N]*ranks[N]),nrow = dims[N])
    Un = qr.Q(qr(T1))
  }
  Dn <- Un%*%S0%*%t(Uj)
  if(d0==0) Dn = array(ttu(Dn,N,1,dims),dims)
  else Dn <- ttu(Dn,N,d0,dims)
 
  return (Dn)
}