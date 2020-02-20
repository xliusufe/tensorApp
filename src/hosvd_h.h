#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <Eigen/SVD>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <list>
#define MIN(a, b) (a<b?a:b)
#define MAX(a, b) (a>b?a:b) 
#define IDEX(a, b) (a>b?1:0) 
using namespace Rcpp;
using namespace Eigen;

//----------------------------------------------------------------**
//***----------------------parameters for penalization---------------**
struct Options
{
	int N;
	double eps;	
	int max_step;
	double eps1;	
	int max_step1;	
	int isfixr;
}opts;

struct Options_pen
{
	int pen; 
	int nlam;
	int dfmax;
	int isPenColumn;
	double lam_max;
	double lam_min;
	double alpha;
	double gamma;
	double eps;
	double eps1;
	double max_step;
	double max_step1;
}opts_pen;

//----------------------------------------------------------------**
//***----------------------Max of a Vector------------------------**
double MaxAbsVector(VectorXd V){
	double x=V[0];
	int n = V.size();
	for(int i=1;i<n;i++)
		x = MAX(x,fabs(V[i]));
	return x;	
}
//----------------------------------------------------------------**
//***----------------------Max of a Vector------------------------**
int MinVectorInt(VectorXi V, int n){
	int x=V[0];
	for(int i=1;i<n;i++)
		x = MIN(x,V[i]);
	return x;	
}
//----------------------------------------------------------------**
//***----------------------cbind----------------------------------**
MatrixXd cbind_rcpp(MatrixXd A, MatrixXd B)
{
	int n = A.rows();
	int p1 = A.cols();
	int p2 = B.cols();
	MatrixXd C = MatrixXd::Constant(n, p1 + p2, 0);
	C.block(0, 0, n, p1) = A;
	C.block(0, p1, n, p2) = B;
	return C;
}
//----------------------------------------------------------------**
//***----------------------rbind----------------------------------**
MatrixXd rbind_rcpp(MatrixXd A, MatrixXd B)
{
	int n1 = A.rows();
	int n2 = B.rows();
	int p = A.cols();
	MatrixXd C = MatrixXd::Constant(n1 + n2, p, 0);
	C.block(0, 0, n1, p) = A;
	C.block(n1, 0, n2, p) = B;
	return C;
}
//----------------------------------------------------------------**
//***--------------------transfer modal 1 to 2 -------------------**
MatrixXd TransferModalUnfoldingsT12(MatrixXd S, VectorXi dim)
{
    int k,d=1, order = dim.size(), r1=dim[0], r2=dim[1];
	for(k=2; k<order; k++) d*=dim[k];
	MatrixXd S1 = S.block(0, 0, r1, r2).transpose();
	for(k=1; k<d; k++) S1 = cbind_rcpp(S1, S.block(0, k*r2, r1, r2).transpose());
	return S1;	
}
//----------------------------------------------------------------**
//***--------------------transfer modal 2 to 1 -------------------**
MatrixXd TransferModalUnfoldingsT21(MatrixXd S, VectorXi dim)
{
    int k,d=1, order = dim.size(), r1=dim[0], r2=dim[1];
	for(k=2; k<order; k++) d*=dim[k];
	MatrixXd S1 = S.block(0, 0, r2, r1).transpose();
	for(k=1; k<d; k++) S1 = cbind_rcpp(S1, S.block(0, k*r1, r2, r1).transpose());
	return S1;	
}
//----------------------------------------------------------------**
//***--------------------transfer modal 1 to d -------------------**
MatrixXd TransferModalUnfoldingsT1d(MatrixXd S, int d, VectorXi dim)
{
	if(d==1) return S;
	if(d==2) return TransferModalUnfoldingsT12(S, dim);
	else{
		d=d-1;
		int i,ii,j,jd,k,d1=1, d2=1, order = dim.size(),r1=dim[0], rd=dim[d];
		for(k=1; k<d; k++) d1*=dim[k];
		for(k=d+1; k<order; k++) d2*=dim[k];		
		MatrixXd S1 ,C = MatrixXd::Constant(r1, rd, 0), C1;
		
		for(ii=0;ii<rd;ii++) C.col(ii) = S.col(d1*ii);
		C1 = C.transpose();
		for(i=1;i<d1;i++){
			for(ii=0;ii<rd;ii++) C.col(ii) = S.col(d1*ii+i);
			C1 = cbind_rcpp(C1,C.transpose());
		}
		S1 = C1;
		for(j=1;j<d2;j++){
			jd = j*d1*rd;
			for(ii=0;ii<rd;ii++) C.col(ii) = S.col(jd+d1*ii);
			C1 = C.transpose();
			for(i=1;i<d1;i++){
				for(ii=0;ii<rd;ii++) C.col(ii) = S.col(jd+d1*ii+i);
				C1 = cbind_rcpp(C1,C.transpose());
			}
			S1 = cbind_rcpp(S1, C1);
			
		}
		return S1;
	}
}
//----------------------------------------------------------------**
//***--------------------transfer modal d to 1 -------------------**
MatrixXd TransferModalUnfoldingsTd1(MatrixXd S, int d, VectorXi dim)
{
	if(d==1) return S;
	if(d==2) return TransferModalUnfoldingsT21(S, dim);
	else{
		d=d-1;
		int i,ii,j,jd,k,d1=1, d2=1, order = dim.size(),r1=dim[0], rd=dim[d];
		for(k=1; k<d; k++) d1*=dim[k];
		for(k=d+1; k<order; k++) d2*=dim[k];		
		MatrixXd S1, C = MatrixXd::Constant(r1, d1*rd, 0), C0;
		for(i=0; i<d1; i++){
			C0 = S.block(0, i*r1, rd, r1).transpose();
			for(ii=0;ii<rd;ii++) C.col(d1*ii+i) = C0.col(ii);				
		}
		S1 = C;
		for(j=1;j<d2;j++){
			jd = j*d1*r1;
			for(i=0; i<d1; i++){
				C0 = S.block(0, jd+i*r1, rd, r1).transpose();
				for(ii=0;ii<rd;ii++) C.col(d1*ii+i) = C0.col(ii);				
			}
			S1 = cbind_rcpp(S1, C);
		}
		return S1;
	}
}
//----------------------------------------------------------------**
//***--------------------transfer modal d1 to d2 -----------------**
// [[Rcpp::export]]
MatrixXd TransferModalUnfoldingsT(MatrixXd T, int d1, int d2, VectorXi dims)
{
	if(dims.size()<3) stop("T must be greater than 3!");
	if(d1==1)
		return TransferModalUnfoldingsT1d(T, d2, dims);
	else if(d2==1) 
		return TransferModalUnfoldingsTd1(T, d1, dims);
	else
		return TransferModalUnfoldingsT1d(
		       TransferModalUnfoldingsTd1(T, d1, dims),
			   d2, dims);
}
//----------------------------------------------------------------**
//***--------------generate a semi-symmetric tensor---------------**
// [[Rcpp::export]]
MatrixXd gtsem0(MatrixXd S, int r1, int r2, VectorXi dims){
	int j,k,i0,flag=1,d=S.size(),count=1,d1=dims[r1-1];
	VectorXi resid;
	MatrixXd S1,S2;
	S1.setZero(d,1);
	resid.setZero(d);
	S.resize(d,1);
	
	for(j=1;j<d;j++){
		S1(j,0) = j;
		resid[j] = j;
	}
	S1.resize(d1,d/d1);
	S2 = TransferModalUnfoldingsT(S1,r1,r2,dims);
	S2.resize(d,1);
	
	while(count<d){
		for(k=count;k<d;k++){
			if(resid[k]>0){
				j = resid[k];
				count = k+1;
				break;
			}
		}
		if(k>=d) break;
		
		i0 = j;
		while(flag){
			if(S2(j,0)==i0)  break;
			else{
				j = S2(j,0);
				S(j,0) = S(i0,0);
				resid[j] = 0;
			}
		}
	}
	S.resize(d1,d/d1);
	return S;
}
//----------------------------------------------------------------**
//***--------------------penalty----------------------------------**
double penalties(double z, double v, double lambda, double alpha, double gamma, int penalty) {
	double beta=0,l1,l2;
	l1 = lambda*alpha; 
	l2 = lambda*(1-alpha);
	if (penalty==1){			  
		if (z > l1) beta = (z-l1)/(v*(1+l2));
		if (z < -l1) beta = (z+l1)/(v*(1+l2));
	}
	if (penalty==2){
		double s = 0;
		if (z > 0) s = 1;
		else if (z < 0) s = -1;
		if (fabs(z) <= l1) beta = 0;
		else if (fabs(z) <= gamma*l1*(1+l2)) beta = s*(fabs(z)-l1)/(v*(1+l2-1/gamma));
		else beta = z/(v*(1+l2));
	}
	if (penalty==3){
		double s = 0;
		if (z > 0) s = 1;
		else if (z < 0) s = -1;
		if (fabs(z) <= l1) beta = 0;
		else if (fabs(z) <= (l1*(1+l2)+l1)) beta = s*(fabs(z)-l1)/(v*(1+l2));
		else if (fabs(z) <= gamma*l1*(1+l2)) beta = s*(fabs(z)-gamma*l1/(gamma-1))/(v*(1-1/(gamma-1)+l2));
		else beta = z/(v*(1+l2));
	}
	return(beta);
}
//----------------------------------------------------------------**
//***----update the jth row of matrix A with penalty--------------**
VectorXd updateAj(VectorXd z, double lambda, double alpha, double gamma, int penalty)
{
	int i,d=z.size();
	for(i=0;i<d;i++)
	z[i] = penalties(z[i], 1, lambda, alpha, gamma, penalty);
	return z;
}
//----------------------------------------------------------------**
//***-----------setup tuning parameters for SCPTPM----------------**
// [[Rcpp::export]]
VectorXd setuplambdaPC(MatrixXd T0, int d0, VectorXi dims, List D0, int nlam, VectorXd setlam)
{
	int N = dims.size(),j,m;
	double lam_max, lam_min, alpha, max_lam, max_tmp,lamj;
	VectorXd lambda, lambda1, Gamma,tmp,Uj;	
	MatrixXd T,T1;
	
	T1 = TransferModalUnfoldingsT(T0,d0,N,dims);
	
	max_tmp=0;
	for(m=0; m<N; m++){	
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
		Uj = T*Gamma;
		lamj = Uj.norm();
		max_tmp = MAX(max_tmp,MaxAbsVector(Uj/lamj));	
	}

	lam_max = setlam[0];
	lam_min = setlam[1];
	alpha = setlam[2];
	max_lam = lam_max * max_tmp / alpha;
	if (lam_min == 0) {
		lambda1.setLinSpaced(nlam - 1, log(max_lam), log(0.0001*max_lam));
		lambda.setLinSpaced(nlam, 0, 0);
		lambda.segment(0, nlam - 1) = lambda1.array().exp().matrix();
	}
	else {
		lambda1.setLinSpaced(nlam, log(max_lam), log(lam_min*max_lam));
		lambda = lambda1.array().exp();
	}

	return lambda;
}