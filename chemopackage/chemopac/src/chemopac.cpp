#include<RcppArmadillo.h>

using namespace Rcpp;


RcppExport SEXP rowSumscpp(SEXP D) {

Rcpp::NumericMatrix Dr(D);                 

int n = Dr.nrow();
int p = Dr.ncol();

arma::mat d(Dr.begin(), n, p, false);    

arma::vec OUT = arma::zeros(n);

for (int i=0; i<n; i++) {
for (int j=0; j<p; j++) {
OUT.at(i) += d.at(i,j);
}}

return as<NumericVector>(wrap(OUT));

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


RcppExport SEXP d2Mhat5cpp(SEXP RESIDUAL, SEXP iaaA, SEXP D, SEXP D2, SEXP R) {

Rcpp::NumericMatrix RESIDUALr(RESIDUAL);                 
Rcpp::NumericMatrix iaaAr(iaaA);                 
Rcpp::NumericMatrix Dr(D);                 
Rcpp::NumericMatrix D2r(D2);                 
Rcpp::NumericVector RRr(R);                 

int n = Dr.nrow();
int na1 = 155;

arma::mat residual(RESIDUALr.begin(), RESIDUALr.nrow(), RESIDUALr.ncol(), false);    
arma::mat iA(iaaAr.begin(), n, n, false);    
arma::mat d(Dr.begin(), n, Dr.ncol(), false);    
arma::mat d2(D2r.begin(), n, 5, false);    
arma::vec r(RRr.begin(), 5, false);    

arma::mat retvec(na1*25,3);
retvec.zeros();
long double V1;
long double V2;
long double V3;
arma::mat B(5,5);
B.eye();

for (int k1=0; k1<5; k1++){
for (int k2=0; k2<5; k2++){

int Fp = na1*(k1*5+k2);
int Lp = na1*(k1*5+k2+1)-1;

for (int i=0; i<n; i++){
int First = na1*i;
int Last = na1*(i+1)-1;
V2 = 0;
for (int j=0; j<n; j++) {
V1 = -d2(j,0)*r(0) - d2(j,1)*r(1) - d2(j,2)*r(2) - d2(j,3)*r(3) - d2(j,4)*r(4);
V1 = exp(V1);
V3 = d(j,k1);
V3 *= d(j,k2);
V3 *= r(k1);
V3 *= r(k2);
V3 *= 4;
V3 -= 2*B(k1,k2)*r(k1);
V1 *= V3;
V1 *= iA(i,j);
V2 += V1;}
retvec(arma::span(Fp,Lp),arma::span::all) += V2*residual(arma::span(First,Last),arma::span::all);}}}
 
return as<NumericMatrix>(wrap(retvec));

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


RcppExport SEXP dMhat5cpp(SEXP RESIDUAL, SEXP iaaA, SEXP D, SEXP D2, SEXP R) {

Rcpp::NumericMatrix RESIDUALr(RESIDUAL);                 
Rcpp::NumericMatrix iaaAr(iaaA);                 
Rcpp::NumericMatrix Dr(D);                 
Rcpp::NumericMatrix D2r(D2);                 
Rcpp::NumericVector RRr(R);                 

int n = Dr.nrow();
int na1 = 155;

arma::mat residual(RESIDUALr.begin(), RESIDUALr.nrow(), RESIDUALr.ncol(), false);    
arma::mat iA(iaaAr.begin(), n, n, false);    
arma::mat d(Dr.begin(), n, Dr.ncol(), false);    
arma::mat d2(D2r.begin(), n, 5, false);    
arma::vec r(RRr.begin(), 5, false);    

arma::mat retvec(na1*5,3);
retvec.zeros();
long double V1;
long double V2;

for (int k=0; k<5; k++) {
int Fp = na1*k;
int Lp = na1*(k+1)-1;
for (int i=0; i<n; i++) {
int First = na1*i;
int Last = na1*(i+1)-1;
V2 = 0;
for (int j=0; j<n; j++) {
V1 = -d2(j,0)*r(0) - d2(j,1)*r(1) - d2(j,2)*r(2) - d2(j,3)*r(3) - d2(j,4)*r(4);
V1 = exp(V1);
V1 *= d(j,k);
V1 *= r(k); 
V1 *= 2;
V1 *= -1;
V1 *= iA(i,j);
V2 += V1;}
retvec(arma::span(Fp,Lp),arma::span::all) += V2*residual(arma::span(First,Last),arma::span::all);}}

return as<NumericMatrix>(wrap(retvec));

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RcppExport SEXP Mhat5cpp(SEXP RESIDUAL, SEXP iaaA, SEXP D2, SEXP R) {

Rcpp::NumericMatrix RESIDUALr(RESIDUAL);                 
Rcpp::NumericMatrix iaaAr(iaaA);                            
Rcpp::NumericMatrix D2r(D2);                 
Rcpp::NumericVector RRr(R);                 

int n = D2r.nrow();
int na1 = 155;

arma::mat residual(RESIDUALr.begin(), RESIDUALr.nrow(), RESIDUALr.ncol(), false);    
arma::mat iA(iaaAr.begin(), n, n, false);    
arma::mat d2(D2r.begin(), n, 5, false);    
arma::vec r(RRr.begin(), 5, false);    

arma::mat retvec(na1,3);
retvec.zeros();
long double V1;
long double V2;

for (int i=0; i<n; i++) {
int First = na1*i;
int Last = na1*(i+1)-1;
V2 = 0;
for (int j=0; j<n; j++) {
V1 = -d2(j,0)*r(0) - d2(j,1)*r(1) - d2(j,2)*r(2) - d2(j,3)*r(3) - d2(j,4)*r(4);
V1 = exp(V1);
V1 *= iA(i,j);
V2 += V1;}
retvec += V2*residual(arma::span(First,Last),arma::span::all);}

return as<NumericMatrix>(wrap(retvec));

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RcppExport SEXP makeaacpp(SEXP THETA, SEXP DESIGN, SEXP R) {
               
Rcpp::NumericMatrix THETAr(THETA);                 
Rcpp::NumericMatrix DESIGNr(DESIGN);
Rcpp::NumericVector Rr(R);

int N = THETAr.nrow();
int P = DESIGNr.ncol();
int M = DESIGNr.nrow();

arma::mat T(THETAr.begin(), N, P, false);       
arma::mat D(DESIGNr.begin(), M, P, false);       
arma::vec r(Rr.begin(), P, false);

arma::mat OUT(N,M);
for (int i=0; i<N; i++) {
for (int j=0; j<M; j++) {
OUT(i,j) = 0;
for(int k=0; k<P; k++) {
OUT(i,j) -= r(k)*(T(i,k)-D(j,k))*(T(i,k)-D(j,k));}}}

OUT = exp(OUT);

return as<NumericMatrix>(wrap(OUT));

}


////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP diagscpp(SEXP X, SEXP Y) {
               
Rcpp::NumericMatrix xr(X);                 
Rcpp::NumericMatrix yr(Y);

int p = xr.nrow(), n = xr.ncol();

arma::mat matx(xr.begin(), p, n, false);       
arma::mat maty(yr.begin(), n, p, false);       

arma::vec D(p);
for (int k=0; k<p; k++) {
D(k)=0;
for (int j=0; j<n; j++) {
D(k) += matx(k,j)*maty(j,k);}}

return as<NumericVector>(wrap(D));

}

////////////////////////////////////////////////////////////////////////////////////////////////////;

RcppExport SEXP solvecpp(SEXP X) {

Rcpp::NumericMatrix xr(X);                 

int r0 = xr.nrow();
int c0 = xr.ncol();

arma::mat x(xr.begin(), r0, c0, false);   

arma::mat y = inv_sympd(x);

return as<NumericMatrix>(wrap(y));

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


RcppExport SEXP odes5cpp(SEXP t, SEXP y, SEXP p) {

Rcpp::NumericVector tr(t);                 
Rcpp::NumericVector yr(y);
Rcpp::NumericVector pr(p);

arma::vec T(tr.begin(), tr.size(), false);       
arma::vec Y(yr.begin(), yr.size(), false);    
arma::vec P(pr.begin(), pr.size(), false);    

double rat = (1/P(7))-(1/P(8));
double y12 = Y(0)*Y(1);

double r1 = 10000*Y(6);
double r2 = P(0)*exp(-(P(3)/P(6))*rat)*Y(0)*Y(2)*Y(3);
double r3 = P(1)*exp(-(P(4)/P(6))*rat)*y12*Y(4);
double r4 = P(2)*exp(-(P(5)/P(6))*rat)*y12;

arma::vec out = arma::zeros(7);

out(0) = 2*r1 - r4;
out(1) = r2 - r3 - r4;
out(2) = r3 - r2;
out(3) = out(2);
out(4) -= out(2);
out(5) -= out(1);
out(5) -= out(3);
out(6) -= r1;

return as<NumericVector>(wrap(out));

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


RcppExport SEXP makeA5cpp(SEXP D1, SEXP D2, SEXP D3, SEXP D4, SEXP D5, SEXP R) {

Rcpp::NumericMatrix d1r(D1);                 
Rcpp::NumericMatrix d2r(D2); 
Rcpp::NumericMatrix d3r(D3); 
Rcpp::NumericMatrix d4r(D4); 
Rcpp::NumericMatrix d5r(D5); 
Rcpp::NumericVector Rr(R); 

int n = d1r.nrow();

arma::mat d1(d1r.begin(), n, n, false);    
arma::mat d2(d2r.begin(), n, n, false);
arma::mat d3(d3r.begin(), n, n, false);
arma::mat d4(d4r.begin(), n, n, false);
arma::mat d5(d5r.begin(), n, n, false);
arma::vec r(Rr.begin(), 5, false); 

arma::mat OUT = arma::zeros(n,n);

for (int i=0; i<n; i++) {
for (int j=0; j<n; j++) {
OUT(i,j)= -r(0)*d1(i,j) - r(1)*d2(i,j) - r(2)*d3(i,j) - r(3)*d4(i,j) - r(4)*d5(i,j);
}}

return as<NumericMatrix>(wrap(OUT));

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RcppExport SEXP KRON2cpp(SEXP aa0, SEXP aa1, SEXP yy) {

Rcpp::NumericMatrix aa0r(aa0);                 
Rcpp::NumericMatrix aa1r(aa1);
Rcpp::NumericVector yr(yy);

int r0 = aa0r.nrow();
int r1 = aa1r.nrow();
int c0 = aa0r.ncol();
int c1 = aa1r.ncol();

arma::mat a0(aa0r.begin(), r0, c0, false);    
arma::mat a1(aa1r.begin(), r1, c1, false);   
arma::vec y(yr.begin(), yr.size(), false); 

int row_offset1;
int col_offset1;
double factor1;
int index;

arma::vec retvec = arma::zeros(r0*r1);

for (int row_idx0=0; row_idx0<r0; ++row_idx0) {

row_offset1 = row_idx0;		
row_offset1 *= r1;	

for (int row_idx1=0; row_idx1<r1; ++row_idx1) {

index = row_offset1 + row_idx1;

for (int col_idx0=0; col_idx0<c0; ++col_idx0 ) {

col_offset1 = col_idx0;		
col_offset1 *= c1;
factor1  = a0.at(row_idx0,col_idx0);

for (int col_idx1=0; col_idx1<c1; ++col_idx1 ) {

retvec.at(index) += factor1 * a1.at(row_idx1,col_idx1) * y.at( col_offset1 + col_idx1 );
				}
			}
		}
	}

return as<NumericVector>(wrap(retvec));

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RcppExport SEXP KRON3cpp(SEXP aa0, SEXP aa1, SEXP aa2, SEXP yy) {

Rcpp::NumericMatrix aa0r(aa0);                 
Rcpp::NumericMatrix aa1r(aa1);
Rcpp::NumericMatrix aa2r(aa2);
Rcpp::NumericVector yr(yy);

int r0 = aa0r.nrow();
int r1 = aa1r.nrow();
int r2 = aa2r.nrow();
int c0 = aa0r.ncol();
int c1 = aa1r.ncol();
int c2 = aa2r.ncol();

arma::mat a0(aa0r.begin(), r0, c0, false);    
arma::mat a1(aa1r.begin(), r1, c1, false);
arma::mat a2(aa2r.begin(), r2, c2, false);   
arma::vec y(yr.begin(), yr.size(), false); 

int row_offset1;
int row_offset2;
int col_offset1;
int col_offset2;
double factor1;
double factor2;
int index;

arma::vec retvec = arma::zeros(r0*r1*r2);

for (int row_idx0=0; row_idx0<r0; ++row_idx0) {

row_offset1 = row_idx0;		
row_offset1    *= r1;

	for (int row_idx1=0; row_idx1<r1; ++row_idx1) {
			
	row_offset2 = row_offset1 + row_idx1;		
	row_offset2 *= r2;	

		for (int row_idx2=0; row_idx2<r2; ++row_idx2) {

		index = row_offset2 + row_idx2;

			for (int col_idx0=0; col_idx0<c0; ++col_idx0 ) {

			col_offset1 = col_idx0;
			col_offset1 *= c1;
			factor1  = a0.at(row_idx0,col_idx0);

				for (int col_idx1=0; col_idx1<c1; ++col_idx1 ) {
						
						col_offset2 = col_offset1 + col_idx1;		
						col_offset2 *= c2;
						factor2  = factor1 * a1.at(row_idx1,col_idx1);	

							for (int col_idx2=0; col_idx2<c2; ++col_idx2 ) {

							retvec.at(index) += factor2 * a2.at(row_idx2,col_idx2) * y.at(col_offset2 + col_idx2);
						}
					}
				}
			}
		}
	}

return as<NumericVector>(wrap(retvec));

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RcppExport SEXP disccpp(SEXP dstar, SEXP psi, SEXP da, SEXP iaa) {
               
Rcpp::NumericMatrix dstarr(dstar);
Rcpp::NumericMatrix psir(psi);    
Rcpp::NumericMatrix dar(da);
Rcpp::NumericMatrix iaar(iaa);

int mm = dar.nrow();
int NN = dstarr.nrow();
int arows = iaar.nrow();
int dcols = dstarr.ncol();

arma::mat DSTAR(dstarr.begin(), NN, dcols, false);       
arma::mat PSI(psir.begin(), NN, 2, false);  
arma::mat DA(dar.begin(), mm, 2, false);  
arma::mat IAA(iaar.begin(), arows, mm, false);       

arma::mat OUT(NN,4);

arma::vec tiA(mm);
arma::vec tt(mm);
double tiAt;
int a;
int b;
int c;

for (int i=0; i<NN; i++) {

OUT(i,0)=0;
OUT(i,1)=0;
OUT(i,2)=0;

for (int p=0; p<mm; p++) {
tt(p) = -PSI(i,0)*DA(p,0) - PSI(i,1)*DA(p,1);}

tt = exp(tt);

tiAt = 0;
for (int j=0; j<mm; j++) {

tiA(j) = 0;
a = (i*mm)+j;
b = mm + j;
c = mm + mm + j;

for (int k=0; k<mm; k++) {
tiA(j) += tt(k)*IAA(a,k);
}

tiAt += tiA(j)*tt(j);
OUT(i,0) += DSTAR(i,j)*tiA(j);
OUT(i,1) += DSTAR(i,b)*tiA(j);
OUT(i,2) += DSTAR(i,c)*tiA(j);}

OUT(i,3) = 1 - tiAt;}

return as<NumericMatrix>(wrap(OUT));

}


////////////////////////////////////////////////////////////////////////////////////////////////////;





















































































