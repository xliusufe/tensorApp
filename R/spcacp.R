
spcacp <- function(Y=NULL,d0=NULL,dims=NULL,nactive=NULL,dr=6,criteria="BIC",penalty="LASSO",D0=NULL,
                   lambda=NULL,nlam=20,lam_min=1e-4,gamma=2,alpha=1,eps=1e-4,max_step=20){
  
  
  if(is.null(Y)) stop("Tensor Y must not be NULL !")
  if(is.null(dims)){
    warning("Dimension dims is NULL !") 
    dims = dim(Y)
  }
  N <- length(dims)

  if(is.null(d0)){
    dd = 1
    Y = matrix(Y,dims[1])
  }
  else dd = d0
  if(dd<1||dd>N) stop("d0 must be in the set {1,2,..., N} !")
  
  if (penalty == "LASSO") pen <- 1
  if (penalty == "MCP")   pen <- 2 
  if (penalty=="SCAD"){    
    gamma <- 3.7
    pen <- 3
  }  
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
  
  
  if(is.null(D0)){
    set.seed(1)
    D0 = list();
    for(j in 1:N){
      U = rnorm(dims[j])
      D0[[j]] = U/sqrt(sum(U^2))
    }
    D0[[N+1]] = 0
  }
  
  if (is.null(lambda)) {
    if (nlam < 1) stop("nlambda must be at least 1")
    setlam = c(1,lam_min,alpha,nlam)
    lambda = setuplambdaPC(Y,dd,dims,D0,nlam,setlam)
  }
  
  opts = list(eps=eps,max_step=max_step,N=N,eps1=eps,max_step1=max_step)
  opts_pen = list(pen=pen,nlam=nlam,lam_max=1,lam_min=lam_min,gamma=gamma,alpha=alpha) 
  
  if(is.null(nactive))  fit <- SCPTPM(Y,dd,dr,dims,D0,lambda,opts,opts_pen)
  else {
    active = setdiff(1:N,nactive)
    fit <- SCPTPM_part(Y,dd,dr,dims,active,D0,lambda,opts,opts_pen)
  }
  prodd = prod(dims)
  df = rss = rep(0,nlam)
  for(l in 0:(nlam-1)){
    resids = ttu(Y,dd,N,dims)
    S = fit[[l*(N+1)+N+1]]
    df0 = 0
    for(k in 1:dr){
      Uj = fit[[l*(N+1)+1]][,k]
      df0 = df0 + sum(abs(Uj)!=0)
      for(j in 2:(N-1)){
        tmp = fit[[l*(N+1)+j]][,k]
        Uj = kronecker(tmp,Uj)
        df0 = df0 + sum(abs(tmp)!=0)
      }
      tmp = fit[[l*(N+1)+N]][,k]
      df0 = df0 + sum(abs(tmp)!=0)
      resids = resids - S[k]*kronecker(tmp,t(Uj))
    }
    df[l+1] = df0
    rss[l+1] = sum(resids^2)
  }

  rss1 =  prodd*log(rss/prodd)
  bic <- switch (criteria,
                 BIC = rss1 + log(prodd)*df,
                 AIC = rss1 + 2*df,
                 GCV = rss*prodd/(prodd-df)^2
  )

  
  selected = which.min(bic)
  lambda_opt = lambda[selected]
  rss_opt = rss[selected]
  df_opt = df[selected]
  Dn=list()
  for(j in 1:N)  Dn[[j]] = fit[[(selected-1)*(N+1)+j]]
  Dn[[N+1]] = fit[[(selected-1)*(N+1)+N+1]]
  
  Tnew = matrix(0,dims[N],prodd/dims[N])
  ranks = Dn[[N+1]]
  for(k in 1:dr){
    Uj = Dn[[1]][,k]
    for(j in 2:(N-1)){
      tmp = Dn[[j]][,k]
      Uj = kronecker(tmp,Uj)
    }
    tmp = Dn[[N]][,k]
    Tnew = Tnew + ranks[k]*kronecker(tmp,t(Uj))
  }
  
  if(is.null(d0))  Tnew = array(TransferModalUnfoldingsT(Tnew,N,1,dims),dim=dims)
  else Tnew = TransferModalUnfoldingsT(Tnew,N,d0,dims)

  return (list(Tnew=Tnew,
               Tn=Dn,
               ranks=ranks,
               rss = rss_opt,
               lambda = lambda,
               selectedID = selected,
               lambda_opt = lambda_opt,
               df = df
               )
  )
}