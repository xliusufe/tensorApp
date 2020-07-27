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
MatrixXd TuckerALS(MatrixXd T1, int d0, VectorXi dims, VectorXi rs, List D0, List optsList){
	//Tucker approximation via alterating least squares famework
	opts.N = as<int>(optsList["N"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	
	int m,j,N=opts.N,step = 0,rprod=1;
	MatrixXd T,Gamma,tmp,svdu, Sn,S0,Dnew;
	for(m=0;m<N-1;m++)  rprod*=rs[m];
	S0.setZero(rs[N-1],rprod);

	while(step < opts.max_step){
		step++;
		for(m=0;m<N;m++){
			T = TransferModalUnfoldingsT(T1,d0,m+1,dims);			
			if(m==0){
				Gamma = as<MatrixXd>(D0[1]);
				for(j=2;j<N;j++){
					tmp = kroneckerProduct(as<MatrixXd>(D0[j]),Gamma);
					Gamma = tmp;
				}						
			}
			else{
				Gamma = as<MatrixXd>(D0[0]);
				for(j=1;j<N;j++){
					if(j!=m){
						tmp = kroneckerProduct(as<MatrixXd>(D0[j]),Gamma);
						Gamma = tmp;
					}
				}
			}				
			JacobiSVD<MatrixXd> svd(T*Gamma, ComputeThinU | ComputeThinV);
			svdu = svd.matrixU().leftCols(rs[m]);
			D0[m] = svdu;		
		}		
		Sn = svdu.transpose()*T*Gamma;
		if((Sn-S0).norm()/(S0.norm()+1)<opts.eps) break;
		S0 = Sn;
	}	
	D0[N] = Sn;
	Dnew = svdu*Sn*Gamma.transpose(); 
	return Dnew;
}

//----------------------------------------------------------------**
//***------------------CP approximation via TPM-------------------**
// [[Rcpp::export]]
MatrixXd CPTPM(MatrixXd T0, int d0, int d, VectorXi dims, List D0, List optsList){

	opts.N = as<int>(optsList["N"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	
	int m,j,k,N=opts.N,step,rprod=1;
	double lamj=0.0,lamj0;
	MatrixXd T,T1,U1,Tnew,Ttemp;
	VectorXd svdu,Gamma,tmp,lambda;
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
		lamj0 = pow(10, 6);
		while(step < opts.max_step){
			step++;
			for(m=0;m<N;m++){
				T = TransferModalUnfoldingsT(T1,N,m+1,dims);
				if(m==0){
					Gamma = as<VectorXd>(D0[1]);
					for(j=2;j<N;j++){
						tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
						Gamma = tmp;
					}						
				}
				else{
					Gamma = as<VectorXd>(D0[0]);
					for(j=1;j<N;j++){
						if(j!=m){
							tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
							Gamma = tmp;
						}
					}
				}				
				svdu = T*Gamma;
				lamj = svdu.norm();
				svdu /= lamj; 
				D0[m] = svdu;		
			}					
			if(fabs(lamj-lamj0)<opts.eps) break;
			lamj0 = lamj;
		}
		lambda[k] = lamj;
		Ttemp = lamj*kroneckerProduct(svdu,Gamma.transpose());  
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
//***------------------CP approximation via TPM-------------------**
// [[Rcpp::export]]
MatrixXd CPTPMorthogon(MatrixXd T0, int d0, int d, VectorXi dims, List D0, List optsList){

	opts.N = as<int>(optsList["N"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	
	int m,j,k,N=opts.N,step,rprod=1;
	double lamj=0.0, lamj0;
	MatrixXd T,T1,U1,Tnew,Ttemp;
	VectorXd svdu,Gamma,tmp,lambda;
	List Dnew(N+1);
	d = MIN(d,MinVectorInt(dims,N));
	
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
		lamj0 = pow(10, 6);
		while(step < opts.max_step){
			step++;				
			for(m=0;m<N;m++){
				T = TransferModalUnfoldingsT(T1,N,m+1,dims);
				if(m==0){
					Gamma = as<VectorXd>(D0[1]);
					for(j=2;j<N;j++){
						tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
						Gamma = tmp;
					}						
				}
				else{
					Gamma = as<VectorXd>(D0[0]);
					for(j=1;j<N;j++){
						if(j!=m){
							tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
							Gamma = tmp;
						}
					}
				}
				if(k==0) svdu = T*Gamma;
				else{
					U1 = (as<MatrixXd>(Dnew[m])).leftCols(k);					
					Ttemp = T*Gamma;
					svdu = Ttemp - U1*(U1.transpose()*Ttemp);				
				}
				lamj = svdu.norm();
				svdu /= lamj; 
				D0[m] = svdu;		
			}					
			if(fabs(lamj-lamj0)<opts.eps) break;
			lamj0 = lamj;
		}			
		lambda[k] = lamj;
		Ttemp = lamj*kroneckerProduct(svdu,Gamma.transpose());  
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
//***------------------CP approximation via TPM-------------------**
// [[Rcpp::export]]
MatrixXd CPTPM_dr(MatrixXd T0, int d0, int d, VectorXi dims, List D0, List optsList){

	opts.N = as<int>(optsList["N"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.eps1 = as<double>(optsList["eps1"]);	
	opts.max_step1 = as<int>(optsList["max_step1"]);	
	
	int m,j,k,N=opts.N,step,rprod=1;
	double lamj=0.0,lamj0;
	MatrixXd T,T1,U1,Tnew,Ttemp;
	VectorXd svdu,Gamma,tmp,lambda;
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
		lamj0 = pow(10, 6);
		while(step < opts.max_step){
			step++;
			for(m=0;m<N;m++){
				T = TransferModalUnfoldingsT(T1,N,m+1,dims);
				if(m==0){
					Gamma = as<VectorXd>(D0[1]);
					for(j=2;j<N;j++){
						tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
						Gamma = tmp;
					}						
				}
				else{
					Gamma = as<VectorXd>(D0[0]);
					for(j=1;j<N;j++){
						if(j!=m){
							tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
							Gamma = tmp;
						}
					}
				}
				svdu = T*Gamma;
				lamj = svdu.norm();
				svdu /= lamj; 
				D0[m] = svdu;		
			}					
			if(fabs(lamj-lamj0)<opts.eps) break;
			lamj0 = lamj;
		}
		lambda[k] = lamj;
		Ttemp = lamj*kroneckerProduct(svdu,Gamma.transpose());  
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
//----------------------------------------------------------------**
//***------------------CP approximation via TPM-------------------**
// [[Rcpp::export]]
MatrixXd CPTPMorthogon_dr(MatrixXd T0, int d0, int d, VectorXi dims, List D0, List optsList){

	opts.N = as<int>(optsList["N"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.eps1 = as<double>(optsList["eps1"]);	
	opts.max_step1 = as<int>(optsList["max_step1"]);	
	
	int m,j,k,N=opts.N,step,rprod=1;
	double lamj=0.0,lamj0;
	MatrixXd T,T1,U1,Tnew,Ttemp;
	VectorXd svdu,Gamma,tmp,lambda;
	List Dnew(N);
	d = MIN(d,MinVectorInt(dims,N));
	
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
		lamj0 = pow(10, 6);
		while(step < opts.max_step){
			step++;
			for(m=0;m<N;m++){
				T = TransferModalUnfoldingsT(T1,N,m+1,dims);
				if(m==0){
					Gamma = as<VectorXd>(D0[1]);
					for(j=2;j<N;j++){
						tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
						Gamma = tmp;
					}						
				}
				else{
					Gamma = as<VectorXd>(D0[0]);
					for(j=1;j<N;j++){
						if(j!=m){
							tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
							Gamma = tmp;
						}
					}
				}
				if(k==0) svdu = T*Gamma;
				else{
					U1 = (as<MatrixXd>(Dnew[m])).leftCols(k);					
					tmp = T*Gamma;
					svdu = tmp - U1*(U1.transpose()*tmp);
				}					
				lamj = svdu.norm();
				svdu /= lamj; 
				D0[m] = svdu;		
			}					
			if(fabs(lamj-lamj0)<opts.eps) break;
			lamj0 = lamj;
		}
		lambda[k] = lamj;
		Ttemp = lamj*kroneckerProduct(svdu,Gamma.transpose());  
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
//----------------------------------------------------------------**
//***------------------CP approximation via TPM-------------------**
// [[Rcpp::export]]
MatrixXd CPALS(MatrixXd T0, int d0, int d, VectorXi dims, List D0, List optsList){

	opts.N = as<int>(optsList["N"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	
	int m,j,N=opts.N,step,rprod=1;
	MatrixXd T,T1,Tnew,svdu,Gamma,tmp,S0;
	
	for(m=0;m<N-1;m++)  rprod*=dims[m];
	S0.setZero(dims[N-1],rprod);
	T1 = TransferModalUnfoldingsT(T0,d0,N,dims);	

	step = 0;
	while(step < opts.max_step){
		step++;
		for(m=0;m<N;m++){
			T = TransferModalUnfoldingsT(T1,d0,m+1,dims).transpose();
			if(m==0){
				Gamma = as<MatrixXd>(D0[1]);
				for(j=2;j<N;j++){
					tmp = KRP(as<MatrixXd>(D0[j]),Gamma);
					Gamma = tmp;
				}						
			}
			else{
				Gamma = as<MatrixXd>(D0[0]);
				for(j=1;j<N;j++){
					if(j!=m){
						tmp = KRP(as<MatrixXd>(D0[j]),Gamma);
						Gamma = tmp;
					}
				}
			}			
			svdu = (Gamma.colPivHouseholderQr().solve(T)).transpose();
			D0[m] = svdu;		
		}		
		Tnew = svdu*Gamma.transpose();
		if((Tnew-S0).norm()/(S0.norm()+1)<opts.eps) break;
		S0 = Tnew;
	}
	return Tnew;
}
//----------------------------------------------------------------**
//***------------------CP approximation via TPM-------------------**
// [[Rcpp::export]]
MatrixXd CPTPMsym2(MatrixXd T0, int d, int k1, int k2, VectorXi dims, List D0, List optsList){

	opts.N = as<int>(optsList["N"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.eps1 = as<double>(optsList["eps1"]);	
	opts.isfixr = as<int>(optsList["isfixr"]);
	
	int m,j,k,N=opts.N,step,rprod=1;
	double lamj=0.0,lamj0;
	MatrixXd T,U1,Tnew,Ttemp,T2,T1=T0;
	VectorXd lambda,svdu,Gamma,tmp;
	VectorXi dims1;
	List Dnew(N+1);
	if(k1>k2){
		m = k1;
		k1 = k2;
		k2 = m;		
	}
	
	for(m=0;m<N-1;m++)  rprod*=dims[m];
	for(m=0;m<N;m++){
		U1.setZero(dims[m],d);
		Dnew[m] = U1;
	}
	lambda.setZero(d);
	Tnew.setZero(dims[N-1],rprod);
				
    for(k=0;k<d;k++){
		step = 0;
		lamj0 = 1000000;
		while(step < opts.max_step){
			step++;			
			for(m=0;m<N;m++){
				if(m!=k1&&m!=k2){
					T = TransferModalUnfoldingsT(T1,N,m+1,dims);
					if(m==0){
						Gamma = as<VectorXd>(D0[1]);
						for(j=2;j<N;j++){
							tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
							Gamma = tmp;
						}						
					}
					else{
						Gamma = as<VectorXd>(D0[0]);
						for(j=1;j<N;j++){
							if(j!=m){
								tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
								Gamma = tmp;
							}
						}
					}										
					svdu = T*Gamma;
					lamj = svdu.norm();
					svdu /= lamj; 
					D0[m] = svdu;	
				}
                				
			}		
			dims1 = dims;
			T2 = TransferModalUnfoldingsT(T1,N,1,dims);
			for(j=0;j<N;j++){
				if(j!=k1&&j!=k2){	
					T = TransferModalUnfoldingsT(T2,1,j+1,dims1);					
					Gamma = as<VectorXd>(D0[j]);
					dims1[j] = 1;
					T2 = TransferModalUnfoldingsT(Gamma.transpose()*T,j+1,1,dims1);											
				}					
			}
			T = TransferModalUnfoldingsT(T2,1,k1+1,dims1);			
			JacobiSVD<MatrixXd> svd(T, ComputeThinU | ComputeThinV);
			svdu = svd.matrixU().col(0);
			tmp = svd.singularValues();	
			lamj = tmp[0];			
			D0[k1] = svdu;	
			D0[k2] = svdu;	
			if(fabs(lamj-lamj0)<opts.eps) break;
			lamj0 = lamj;
		}// End while				
		Gamma = as<VectorXd>(D0[0]);
		for(m=1;m<N-1;m++){
			tmp = kroneckerProduct(as<VectorXd>(D0[m]),Gamma);
			Gamma = tmp;
		}
		svdu = as<VectorXd>(D0[N-1]);
		lamj = svdu.transpose()*T1*Gamma;
		lambda[k] = lamj;
		Ttemp = lamj*kroneckerProduct(svdu,Gamma.transpose());  
		T1 -= Ttemp;			
		Tnew += Ttemp;
        for(m=0;m<N;m++){
			U1 = as<MatrixXd>(Dnew[m]);
			U1.col(k) = as<VectorXd>(D0[m]);
			Dnew[m] = U1;
		}
		if(!opts.isfixr) if(T1.norm()<opts.eps1) break;	
        	
	}// End for
	for(m=0;m<N;m++){
		U1 = as<MatrixXd>(Dnew[m]);
		D0[m] = U1;
	}	
    if(opts.isfixr) D0[N] = lambda;	
	else{
		if(k>0) D0[N] = lambda.head(k-1);
	}
	return Tnew;
}
//----------------------------------------------------------------**
//***------------------CP approximation via TPM-------------------**
// [[Rcpp::export]]
MatrixXd CPTPMsym2Orth(MatrixXd T0, int d, int k1, int k2, VectorXi dims, List D0, List optsList, List optsList_pen){

	opts.N = as<int>(optsList["N"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	opts.eps1 = as<double>(optsList["eps1"]);	
	opts.isfixr = as<int>(optsList["isfixr"]);
	
	int m,j,k,N=opts.N,step,rprod=1;
	double lamj=0.0,lamj0;
	MatrixXd T,U1,Tnew,Ttemp,T2,T1=T0;
	VectorXd lambda,svdu,Gamma,tmp;
	VectorXi dims1;
	List Dnew(N+1);
	if(k1>k2){
		m = k1;
		k1 = k2;
		k2 = m;		
	}
	
	for(m=0;m<N-1;m++)  rprod*=dims[m];
	for(m=0;m<N;m++){
		U1.setZero(dims[m],d);
		Dnew[m] = U1;
	}
	lambda.setZero(d);
	Tnew.setZero(dims[N-1],rprod);
				
    for(k=0;k<d;k++){
		step = 0;
		lamj0 = 1000000;
		while(step < opts.max_step){
			step++;			
			for(m=0;m<N;m++){
				if(m!=k1&&m!=k2){
					T = TransferModalUnfoldingsT(T1,N,m+1,dims);
					if(m==0){
						Gamma = as<VectorXd>(D0[1]);
						for(j=2;j<N;j++){
							tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
							Gamma = tmp;
						}						
					}
					else{
						Gamma = as<VectorXd>(D0[0]);
						for(j=1;j<N;j++){
							if(j!=m){
								tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
								Gamma = tmp;
							}
						}
					}										
					svdu = T*Gamma;
					lamj = svdu.norm();
					svdu /= lamj; 
					D0[m] = svdu;	
				}
                				
			}		

			dims1 = dims;
			T2 = TransferModalUnfoldingsT(T1,N,1,dims);
			for(j=0;j<N;j++){
				if(j!=k1&&j!=k2){	
					T = TransferModalUnfoldingsT(T2,1,j+1,dims1);					
					Gamma = as<VectorXd>(D0[j]);
					dims1[j] = 1;
					T2 = TransferModalUnfoldingsT(Gamma.transpose()*T,j+1,1,dims1);											
				}					
			}
			T = TransferModalUnfoldingsT(T2,1,k1+1,dims1);			
			
			if(k==0){
				JacobiSVD<MatrixXd> svd(T, ComputeThinU | ComputeThinV);
				svdu = svd.matrixU().col(0);
				tmp = svd.singularValues();	
				lamj = tmp[0];
			}
			else{
				U1 = (as<MatrixXd>(Dnew[k1])).leftCols(k);					
				T2 = T - U1*(U1.transpose()*T);
				Ttemp = T2 - (T2*U1)*U1.transpose();
				JacobiSVD<MatrixXd> svd(Ttemp, ComputeThinU | ComputeThinV);
				svdu = svd.matrixU().col(0);
				tmp = svd.singularValues();	
				lamj = tmp[0];
			}
			D0[k1] = svdu;	
			D0[k2] = svdu;	
			if(fabs(lamj-lamj0)<opts.eps) break;
			lamj0 = lamj;
		} // End while		
		Gamma = as<VectorXd>(D0[0]);
		for(m=1;m<N-1;m++){
			tmp = kroneckerProduct(as<VectorXd>(D0[m]),Gamma);
			Gamma = tmp;
		}
		svdu = as<VectorXd>(D0[N-1]);
		lamj = svdu.transpose()*T1*Gamma;
		lambda[k] = lamj;
		Ttemp = lamj*kroneckerProduct(svdu,Gamma.transpose());  
		T1 -= Ttemp;			
		Tnew += Ttemp;
        for(m=0;m<N;m++){
			U1 = as<MatrixXd>(Dnew[m]);
			U1.col(k) = as<VectorXd>(D0[m]);
			Dnew[m] = U1;
		}
		if(!opts.isfixr) if(T1.norm()<opts.eps1) break;	
        	
	}// End for	
	for(m=0;m<N;m++){
		U1 = as<MatrixXd>(Dnew[m]);
		D0[m] = U1;
	}	
	if(opts.isfixr) D0[N] = lambda;	
	else{
		if(k>0) D0[N] = lambda.head(k-1);
	}
	return Tnew;
}
//----------------------------------------------------------------**
//***------------------Tucker approximation via ALS---------------**
// [[Rcpp::export]]
MatrixXd TuckerALSsym2(MatrixXd T0, int k1, int k2, VectorXi dims, VectorXi rs, List D0, List optsList){
	//Tucker approximation via alterating least squares for semi-symmetric tensor
	opts.N = as<int>(optsList["N"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	
	int m,j,N=opts.N,step = 0,rprod=1;
	MatrixXd T,Gamma,tmp,svdu, Sn,S0,Dnew,T1=T0,T2;
	VectorXi dims1;
	for(m=0;m<N;m++)  if(m!=k2) rprod*=rs[m];
	S0.setZero(rs[k2],rprod);

	while(step < opts.max_step){
		step++;
		for(m=0;m<N;m++){
			if(m!=k1&&m!=k2){
				T = TransferModalUnfoldingsT(T1,N,m+1,dims);			
				if(m==0){
					Gamma = as<MatrixXd>(D0[1]);
					for(j=2;j<N;j++){
						tmp = kroneckerProduct(as<MatrixXd>(D0[j]),Gamma);
						Gamma = tmp;
					}						
				}
				else{
					Gamma = as<MatrixXd>(D0[0]);
					for(j=1;j<N;j++){
						if(j!=m){
							tmp = kroneckerProduct(as<MatrixXd>(D0[j]),Gamma);
							Gamma = tmp;
						}
					}
				}				
				JacobiSVD<MatrixXd> svd(T*Gamma, ComputeThinU | ComputeThinV);
				svdu = svd.matrixU().leftCols(rs[m]);
				D0[m] = svdu;	
			}
		} 		
		dims1 = dims;
		T2 = TransferModalUnfoldingsT(T1,N,1,dims);
		for(j=0;j<N;j++){
			if(j!=k1&&j!=k2){	
				T = TransferModalUnfoldingsT(T2,1,j+1,dims1);					
				Gamma = as<MatrixXd>(D0[j]);
				dims1[j] = rs[j];
				T2 = TransferModalUnfoldingsT(Gamma.transpose()*T,j+1,1,dims1);								
			}					
		}
		T = TransferModalUnfoldingsT(T2,1,k1+1,dims1);
		JacobiSVD<MatrixXd> svd(T, ComputeThinU | ComputeThinV);
		svdu = svd.matrixU().leftCols(rs[k1]);		
		D0[k1] = svdu;	
		D0[k2] = svdu;	
	
		dims1[k1] = rs[k1];		
		T2 = TransferModalUnfoldingsT(svdu.transpose()*T,k1+1,k2+1,dims1);		
		Sn = svdu.transpose()*T2; 	
		
		if((Sn-S0).norm()/(S0.norm()+1)<opts.eps) break;
		S0 = Sn;
	}
	S0 = TransferModalUnfoldingsT(Sn,k2+1,N,rs);
	D0[N] = S0;
	
	Gamma = as<MatrixXd>(D0[0]);
	for(j=1;j<N-1;j++){
		tmp = kroneckerProduct(as<MatrixXd>(D0[j]),Gamma);
		Gamma = tmp;
	}
	svdu = as<MatrixXd>(D0[N-1]);
	Dnew = svdu*S0*Gamma.transpose(); 
	return Dnew;
}

//-----------------------------------------------------------------**
//***-------------sparseCP approximation via TPM-------------------**
// [[Rcpp::export]]
List SCPTPM(MatrixXd T0, int d0, int d, VectorXi dims, List D1, VectorXd lambda, List optsList, List optsList_pen){
	opts.N = as<int>(optsList["N"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	
	opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
	opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
	opts_pen.gamma = as<double>(optsList_pen["gamma"]);
	opts_pen.alpha = as<double>(optsList_pen["alpha"]);
	opts_pen.pen = as<int>(optsList_pen["pen"]);
	opts_pen.nlam = lambda.size();	
	
	int nlam=opts_pen.nlam, penalty = opts_pen.pen;
	double alpha = opts_pen.alpha, eps = opts.eps, gamma = opts_pen.gamma;
	int l,m,j,k,N=opts.N,step;
	double lamj=0.0,lamj0,lambda1;
	MatrixXd T,T1,U1,Ttemp;
	VectorXd svdu,Gamma,tmp,S,svdu1;
	List Dnew(N+1), Dout(nlam*(N+1)),D0(N+1);
	
	for(m=0;m<N;m++){
		U1.setZero(dims[m],d*N);
		Dnew[m] = U1;
	}
	
	S.setZero(d);

    for(l=0;l<nlam;l++){		
		lambda1 = lambda[l];
		T1 = TransferModalUnfoldingsT(T0,d0,N,dims);
		for(m=0;m<N;m++) D0[m] = as<VectorXd>(D1[m]);		
		for(k=0;k<d;k++){				
			step = 0;
			lamj0 = 100000.0;
			while(step < opts.max_step){
				step++;
				for(m=0;m<N;m++){
					T = TransferModalUnfoldingsT(T1,N,m+1,dims);
					if(m==0){
						Gamma = as<VectorXd>(D0[1]);
						for(j=2;j<N;j++){
							tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
							Gamma = tmp;
						}						
					}
					else{
						Gamma = as<VectorXd>(D0[0]);
						for(j=1;j<N;j++){
							if(j!=m){
								tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
								Gamma = tmp;
							}
						}
					}			
					svdu1 = T*Gamma;
					lamj = svdu1.norm();
					svdu1 /= lamj; 
					svdu = updateAj(svdu1, lambda1, alpha, gamma, penalty);
					lamj = svdu.norm();
					if(lamj>1e-15) svdu /= lamj;
					D0[m] = svdu;						
				}				
				if(fabs(lamj-lamj0)<eps) break;
				lamj0 = lamj;
			}//end while
			Gamma = as<VectorXd>(D0[0]);
			for(m=1;m<N-1;m++){
				tmp = kroneckerProduct(as<VectorXd>(D0[m]),Gamma);
				Gamma = tmp;
			}
			svdu = as<VectorXd>(D0[N-1]);
			lamj = svdu.transpose()*T1*Gamma;			
			S[k] = lamj;
			Ttemp = lamj*kroneckerProduct(svdu,Gamma.transpose());  
			T1 -= Ttemp;
			for(m=0;m<N;m++){
				U1 = as<MatrixXd>(Dnew[m]);
				U1.col(k) = as<VectorXd>(D0[m]);
				Dnew[m] = U1;
			}
						
		}// end for(k=0;k<d;k++) 
		for(m=0;m<N;m++){
			U1 = as<MatrixXd>(Dnew[m]);
			Dout[l*(N+1)+m] = U1;
		}			
		Dout[l*(N+1)+N] = S;	
	}// end for(l=0;l<nlam;l++)
	return Dout;
}
//-----------------------------------------------------------------**
//***-------------sparseCP approximation via TPM-------------------**
// [[Rcpp::export]]
List SCPTPM_part(MatrixXd T0, int d0, int d, VectorXi dims, VectorXi actives, List D1, VectorXd lambda, List optsList, List optsList_pen){
	opts.N = as<int>(optsList["N"]);
	opts.eps = as<double>(optsList["eps"]);	
	opts.max_step = as<int>(optsList["max_step"]);
	
	opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
	opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
	opts_pen.gamma = as<double>(optsList_pen["gamma"]);
	opts_pen.alpha = as<double>(optsList_pen["alpha"]);
	opts_pen.pen = as<int>(optsList_pen["pen"]);
	opts_pen.nlam = lambda.size();	
	
	int nlam=opts_pen.nlam, penalty = opts_pen.pen;
	double alpha = opts_pen.alpha, eps = opts.eps, gamma = opts_pen.gamma;
	int l,m,i,j,k,N=opts.N,step;
	double lamj=0.0,lamj0,lambda1, nas = actives.size(), nnas = N-nas;
	MatrixXd T,T1,U1,Ttemp;
	VectorXd svdu,Gamma,tmp,S,svdu1;
	VectorXi nonactives;
	nonactives.setZero(nnas);
	List Dnew(N+1), Dout(nlam*(N+1)),D0(N+1);
	
	for(m=0;m<N;m++){
		U1.setZero(dims[m],d*N);
		Dnew[m] = U1;
	}
	
	S.setZero(d);

    for(l=0;l<nlam;l++){		
		lambda1 = lambda[l];
		T1 = TransferModalUnfoldingsT(T0,d0,N,dims);
		for(m=0;m<N;m++) D0[m] = as<VectorXd>(D1[m]);		
		for(k=0;k<d;k++){				
			step = 0;
			lamj0 = 100000.0;
			while(step < opts.max_step){
				step++;
				for(i=0;i<nnas;i++){
					m = nonactives[i];
					T = TransferModalUnfoldingsT(T1,N,m+1,dims);
					if(m==0){
						Gamma = as<VectorXd>(D0[1]);
						for(j=2;j<N;j++){
							tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
							Gamma = tmp;
						}						
					}
					else{
						Gamma = as<VectorXd>(D0[0]);
						for(j=1;j<N;j++){
							if(j!=m){
								tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
								Gamma = tmp;
							}
						}
					}					
                    svdu = T*Gamma;
					lamj = svdu.norm();
					if(lamj>1e-15) svdu /= lamj; 
					D0[m] = svdu;		
				}
				
				for(i=0;i<nas;i++){
					m = actives[i];
					T = TransferModalUnfoldingsT(T1,N,m+1,dims);
					if(m==0){
						Gamma = as<VectorXd>(D0[1]);
						for(j=2;j<N;j++){
							tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
							Gamma = tmp;
						}						
					}
					else{
						Gamma = as<VectorXd>(D0[0]);
						for(j=1;j<N;j++){
							if(j!=m){
								tmp = kroneckerProduct(as<VectorXd>(D0[j]),Gamma);
								Gamma = tmp;
							}
						}
					}			
					svdu1 = T*Gamma;
					lamj = svdu1.norm();
					svdu1 /= lamj; 
					svdu = updateAj(svdu1, lambda1, alpha, gamma, penalty);
					lamj = svdu.norm();
					if(lamj>1e-10) svdu /= lamj;
					D0[m] = svdu;						
				}				
				if(fabs(lamj-lamj0)<eps) break;
				lamj0 = lamj;
			}//end while
			Gamma = as<VectorXd>(D0[0]);
			for(m=1;m<N-1;m++){
				tmp = kroneckerProduct(as<VectorXd>(D0[m]),Gamma);
				Gamma = tmp;
			}
			svdu = as<VectorXd>(D0[N-1]);
			lamj = svdu.transpose()*T1*Gamma;			
			S[k] = lamj;
			Ttemp = lamj*kroneckerProduct(svdu,Gamma.transpose());  
			T1 -= Ttemp;
			for(m=0;m<N;m++){
				U1 = as<MatrixXd>(Dnew[m]);
				U1.col(k) = as<VectorXd>(D0[m]);
				Dnew[m] = U1;
			}						
		}// end for(k=0;k<d;k++) 
		for(m=0;m<N;m++){
			U1 = as<MatrixXd>(Dnew[m]);
			Dout[l*(N+1)+m] = U1;
		}				
		Dout[l*(N+1)+N] = S;	
	}// end for(l=0;l<nlam;l++)
	return Dout;
}