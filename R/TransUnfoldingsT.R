
TransUnfoldingsT <- function(S,d1=NULL,d2=0,dims=NULL){
  
  if(is.null(S)) stop("Tensor T must not be NULL !")
  if(is.null(dims)){
    warning("Dimension dims is NULL !") 
    dims = dim(S)
  }
  N <- length(dims)
  if(N<3){
    warning("The input S is not a tensor with dimension greater than 2! It is a matrix.")
    return(S)
  }
  if(is.null(d1)){
    dd = 1
    S = matrix(S,dims[1])
  }
  else dd = d1
  if(dd<1||dd>N) stop("d0 must be in the set {1,2,..., N} !")
  
  if(d2==0){
    Td2 = array(TransferModalUnfoldingsT(S,dd,1,dims),dim=dims)
  }
  else if(d2==1) Td2 = S
  else  Td2 = TransferModalUnfoldingsT(S,dd,d2,dims)
  
  return (Td2)
}