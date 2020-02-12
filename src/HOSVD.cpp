//[[Rcpp::depends(RcppEigen)]]
#include "hosvd_h.h"

//----------------------------------------------------------------**
//***--------------------Estimation without penalty---------------**
// [[Rcpp::export]]
MatrixXd KRP(MatrixXd A, MatrixXd B){//Khatri-Rao Product
	int i,m=A.rows(), n=B.rows(), p=A.cols();
	if(p!=B.cols()) stop("the number of columns both A and B must be same!");
	MatrixXd out;
	out.setZero(m*n,p);
	for(i=0;i<p;i++){
		out.col(i) = kroneckerProduct(A.col(i), B.col(i));
	}
	return out;
}
//----------------------------------------------------------------**
//***------------------Tucker approximation via ALS---------------**
// [[Rcpp::export]]
MatrixXd TuckerALS(MatrixXd T1, int d0, VectorXi dims, VectorXi rs, List Dn, List optsList){
	//Tucker approximation via alterating least squares famewor
	opts.N = as<int>(optsList["N"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	
	int m,j,N=opts.N,step = 0;
	MatrixXd T,Gamma,Gamma1,Gamma2,tmp,svdu,U0, Sn,Dnew;
	VectorXi convergence;
	convergence.setOnes(N);

	while(step < opts.max_step){
		step++;
		for(m=0;m<N;m++){
			T = TransferModalUnfoldingsT(T1,d0,m+1,dims);
			U0 = as<MatrixXd>(Dn[m]);
			if(m!=0) Gamma1 = as<MatrixXd>(Dn[0]);
			for(j=1;j<m;j++){
				tmp = kroneckerProduct(as<MatrixXd>(Dn[j]),Gamma1);
				Gamma1 = tmp;
			}				
			if(m!=N-1)Gamma2 = as<MatrixXd>(Dn[m+1]);
			for(j=m+2;j<N;j++){
				tmp = kroneckerProduct(as<MatrixXd>(Dn[j]),Gamma2);
				Gamma2 = tmp;
			}				
			if(m==0) Gamma = Gamma2;
			else if (m==N-1) Gamma = Gamma1;
			else  Gamma = kroneckerProduct(Gamma2,Gamma1);	
			JacobiSVD<MatrixXd> svd(T*Gamma, ComputeThinU | ComputeThinV);
			svdu = svd.matrixU().leftCols(rs[m]);
			if((svdu-U0).norm()/U0.norm()<opts.eps) convergence[m] = 0;
			Dn[m] = svdu;		
		}		
		if(convergence.sum()==0) break;
	}	
	Sn = svdu.transpose()*T*Gamma;
	Dn[N] = Sn;
	Dnew = svdu*Sn*Gamma.transpose(); 
	return Dnew;
}

//----------------------------------------------------------------**
//***------------------CP approximation via ALS-------------------**
// [[Rcpp::export]]
MatrixXd CPALS(MatrixXd T0, int d0, int d, VectorXi dims, List D0, List optsList){
	/*
	CANDECOMP/PARAFAC Decomposition approximation via alterating least squares famewor
	References:
	Allen, G., 2012. 
	Sparse higher-order principal components analysis, 
	in: International Conference on Artificial Intelligence and Statistics, pp. 27-36.
	
	Zhengwu Zhanga, Genevera I. Allenb,c, Hongtu Zhud, David Dunson (2018).
	Tensor network factorizations: Relationships between brain structural connectomes and traits.
	Neuroimage
	*/
	opts.N = as<int>(optsList["N"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	
	int m,j,k,N=opts.N,step,rprod=1;
	MatrixXd T,T1,U1,Tnew,Ttemp;
	VectorXd svdu, U0,Gamma,Gamma1,Gamma2,tmp,lambda;
	VectorXi convergence;
	List Dnew(N+1);
	
	for(m=0;m<N-1;m++)  rprod*=dims[m];
	for(m=0;m<N;m++){
		U1.setZero(dims[m],d);
		Dnew[m] = U1;
	}
	T1 = TransferModalUnfoldingsT(T0,d0,N,dims);
	lambda.setZero(d);
	Tnew.setZero(dims[N-1],rprod);

    for(k=0;k<d;k++){
		step = 0;
		convergence.setOnes(N);
		while(step < opts.max_step){
			step++;
			for(m=0;m<N;m++){
				T = TransferModalUnfoldingsT(T1,N,m+1,dims);
				U0 = as<VectorXd>(D0[m]);
				if(m!=0) Gamma1 = as<VectorXd>(D0[0]);
				for(j=1;j<m;j++){
					tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma1);
					Gamma1 = tmp;
				}				
				if(m!=N-1)Gamma2 = as<VectorXd>(D0[m+1]);
				for(j=m+2;j<N;j++){
					tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma2);
					Gamma2 = tmp;
				}				
				if(m==0) Gamma = Gamma2;
				else if (m==N-1) Gamma = Gamma1;
				else  Gamma = kroneckerProduct(Gamma2,Gamma1);	
				svdu = T*Gamma;
				svdu /=svdu.norm(); 
				if((svdu-U0).norm()/U0.norm()<opts.eps) convergence[m] = 0;
				D0[m] = svdu;		
			}					
			if(convergence.sum()==0) break;
		}
		lambda[k] = svdu.transpose()*T*Gamma;
		Ttemp = lambda[k]*kroneckerProduct(svdu,Gamma.transpose());  
		T1 -= Ttemp;			
		Tnew += Ttemp;
        for(m=0;m<N;m++){
			U1 = as<MatrixXd>(Dnew[m]);
			U1.col(k) = as<VectorXd>(D0[m]);
			Dnew[m] = U1;
		}		
	}
	for(m=0;m<N;m++){
		U1 = as<MatrixXd>(Dnew[m]);
		D0[m] = U1;
	}	
    D0[N] = lambda;	
	return Tnew;
}
//----------------------------------------------------------------**
//***------------------CP approximation via ALS-------------------**
// [[Rcpp::export]]
MatrixXd CPALS_dr(MatrixXd T0, int d0, int d, VectorXi dims, List D0, List optsList){
	/*
	CANDECOMP/PARAFAC Decomposition approximation via alterating least squares famewor
	References:
	Allen, G., 2012. 
	Sparse higher-order principal components analysis, 
	in: International Conference on Artificial Intelligence and Statistics, pp. 27-36.
	
	Zhengwu Zhanga, Genevera I. Allenb,c, Hongtu Zhud, David Dunson (2018).
	Tensor network factorizations: Relationships between brain structural connectomes and traits.
	Neuroimage
	*/
	opts.N = as<int>(optsList["N"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.eps1 = as<double>(optsList["eps1"]);	
	opts.max_step1 = as<int>(optsList["max_step1"]);	
	
	int m,j,k,N=opts.N,step,rprod=1;
	MatrixXd T,T1,U1,Tnew,Ttemp;
	VectorXd svdu, U0,Gamma,Gamma1,Gamma2,tmp,lambda;
	VectorXi convergence;
	List Dnew(N);
	for(m=0;m<N;m++){
		U1.setZero(dims[m],d);
		Dnew[m] = U1;
	}
	for(m=0;m<N-1;m++)  rprod*=dims[m];	
	
	T1 = TransferModalUnfoldingsT(T0,d0,N,dims);
	lambda.setZero(d);
	Tnew.setZero(dims[N-1],rprod);
	
    for(k=0;k<d;k++){
		step = 0;
		convergence.setOnes(N);
		while(step < opts.max_step){
			step++;
			for(m=0;m<N;m++){
				T = TransferModalUnfoldingsT(T1,N,m+1,dims);
				U0 = as<VectorXd>(D0[m]);
				if(m!=0) Gamma1 = as<VectorXd>(D0[0]);
				for(j=1;j<m;j++){
					tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma1);
					Gamma1 = tmp;
				}				
				if(m!=N-1)Gamma2 = as<VectorXd>(D0[m+1]);
				for(j=m+2;j<N;j++){
					tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma2);
					Gamma2 = tmp;
				}				
				if(m==0) Gamma = Gamma2;
				else if (m==N-1) Gamma = Gamma1;
				else  Gamma = kroneckerProduct(Gamma2,Gamma1);	
				svdu = T*Gamma;
				svdu /=svdu.norm(); 
				if((svdu-U0).norm()/U0.norm()<opts.eps) convergence[m] = 0;
				D0[m] = svdu;		
			}					
			if(convergence.sum()==0) break;
		}
		lambda[k] = svdu.transpose()*T*Gamma;
		Ttemp = lambda[k]*kroneckerProduct(svdu,Gamma.transpose());  
		T1 -= Ttemp;			
		Tnew += Ttemp;		
		
        for(m=0;m<N;m++){
			U1 = as<MatrixXd>(Dnew[m]);
			U1.col(k) = as<VectorXd>(D0[m]);
			Dnew[m] = U1;
		}
		if(T1.norm()<opts.eps1) break;	
	}
	for(m=0;m<N;m++){
		U1 = as<MatrixXd>(Dnew[m]);
		D0[m] = U1;
	}	
	if(k>0)	D0[N] = lambda.head(k-1);
	else D0[N] = lambda.head(0);
	return Tnew;
}