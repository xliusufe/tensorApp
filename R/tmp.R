
tmp <- function(S=NULL,M=NULL,d0=NULL,dims=NULL){
  
  if(is.null(S)) stop("Tensor Y must not be NULL !")
  if(is.null(M)) stop("Matrix M must not be NULL !")
  if(is.null(d0)) stop("Mode d must be given !")
  if(is.null(dims))  dims = dim(S)
  if(dims[d0]!=ncol(M)) stop(paste("The ",as.character(d0),"s dimension must be equal to the number of columns of M !",sep = ""))
  
  N <- length(dims)
  if(d0<1||d0>N) stop("d must be in the set {1,2,..., N} !")
  if(N<3){
    warning("The input S is not a tensor with dimension greater than 2! This is an ordinary multiplication of matrices.")
    return (S%*%M)
  } 
  S = matrix(S,dims[1])
  Tnew = M%*%TransferModalUnfoldingsT(S,1,d0,dims)
  Tnew = array(TransferModalUnfoldingsT(Tnew,d0,1,dims),dims)
  return (Tnew)
}