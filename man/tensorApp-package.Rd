\name{tensorApp-package}
\alias{tensorApp-package}
\alias{tensorApp}
\docType{package}
\title{
  High-order SVD approximation of a tensor \code{Y} by Tucker or CP decomposition and selection of ranks
}
\description{
  High-order SVD approximation of a tensor \code{Y} by Tucker decomposition or CANDECOMP/PARAFAC (CP) decomposition and selection of ranks. 
  Alternating Least Squares algorithm is applied to Tuchker decomposition, and both Alternating Least Squares algorithm or 
  Tensor Power Method are applied to CP decomposition. This package provides several generator functions, which generate low-rank tensor or low-rank semi-symmetric tensor. 
}
\details{
  High-order SVD approximation of a tensor \code{Y} by Tucker decomposition or CANDECOMP/PARAFAC (CP) decomposition and selection of ranks.
}
\author{
Xu Liu

Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
}
\references{
Allen, G., (2012). Sparse higher-order principal components analysis, in: International Conference on Artificial Intelligence
and Statistics, pp. 27-36.

De Lathauwer, L., De Moor, B., Vandewalle, J., (2000). A multilinear singular value decomposition. SIAM Journal on
Matrix Analysis and Applications, \eqn{{\bf 21}}, 1253-1278.

De Lathauwer, L., De Moor, B., Vandewalle, J., (2000). On the best rank-1 and rank-\eqn{(r_1, r_2,\cdots, r_n)} approximation of
higher-order tensors. SIAM Journal on Matrix Analysis and Applications, \eqn{{\bf 21}}, 1324-1342.

Kolda, T.G., (2001). Orthogonal tensor decompositions. SIAM Journal on Matrix Analysis and Applications, \eqn{{\bf 23}}, 243-255.

Kolda, T., Bader, B., (2009). Tensor decompositions and applications. SIAM Review, \eqn{{\bf 51}}, 455-500.
}
\keyword{CP decomposition; HOSVD; Tucker decomposition}
