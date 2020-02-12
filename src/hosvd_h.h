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
}opts;

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
