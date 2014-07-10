/* @file test.cpp */

#include "IOFiles.hpp"
#include "Gnuplot.hpp"
#include "Lapack.hpp"

#include <stdlib.h>
#include <time.h>

std::complex<double> projection(Matrix<double> const& O, Matrix<std::complex<double> > const& base, unsigned int bra, unsigned int ket);
Matrix<double> create_to_check_constructor();

int main(){
	std::complex<double> a(-6e-16,23);
	std::complex<double> b(5e-16,23);
	std::cout<<(are_equal(a,b,1e-14,1e-13)?"a==b":"a!=b")<<std::endl;
	/*{operateurs*/
	//unsigned int N1_row(3);
	//unsigned int N1_col(2);
	////unsigned int N2_row(2);
	////unsigned int N2_col(6);
	//Matrix<double> m1(N1_row,N1_col);
	////Matrix<double> m2(N2_row,N2_col,2);
	////
	//m1(0,0) = 1;
	//m1(0,1) = 2;
	//m1(1,0) = 3;
	//m1(1,1) = 4;
	//m1(2,0) = 5;
	//m1(2,1) = 6;
	////Matrix<double> m3(m1);
	////
	////std::cout<<"m1"<<std::endl;
	////std::cout<<m1<<std::endl;
	////std::cout<<"m2"<<std::endl;
	////std::cout<<m2<<std::endl;
	////std::cout<<"m3"<<std::endl;
	////std::cout<<m3<<std::endl;
	////m1 = m1*m2;
	////std::cout<<"m1=m1*m2"<<std::endl;
	////std::cout<<m1<<std::endl;
	////m2=m2.transpose();
	////std::cout<<"m2=m2.tanspose"<<std::endl;
	////std::cout<<m2<<std::endl;
	////std::cout<<"m1*m2"<<std::endl;
	////std::cout<<m1*m2<<std::endl;
	////std::cout<<"m3-(m1-m2)"<<std::endl;
	////std::cout<<m3-(m1*m2)<<std::endl;
	////
	////Matrix<double> m4(1,3,3.1);
	////Matrix<double> m5(3,1,11);
	////
	////std::cout<<"m4*m5"<<std::endl;
	////std::cout<<m5*m4<<std::endl;
	////m4 *= m1;
	////std::cout<<m4<<std::endl;
	////m1.print_mathematica();
	////std::cout<<m1<<std::endl;
	//Vector<double> v(m1.col());
	//for(unsigned int i(0);i<v.size();i++){v(i) = i+1;} 
	//std::cout<<m1<<std::endl;
	//std::cout<<std::endl;
	//std::cout<<v<<std::endl;
	//std::cout<<std::endl;
	//std::cout<<m1*v<<std::endl;
	/*}*/
	///*{eigenvalue*/
	///*{real symmetric matrix*/
	//srand(time(NULL));
	//unsigned int N_site(11);
	//Matrix<double> S(N_site,N_site,0.0);
	//S(0,N_site-1)=1.0;
	//for(unsigned int i(0);i<S.row();i++){
		//for(unsigned int j(i);j<S.col();j++){
			//S(i,j) = (1.0*rand())/RAND_MAX;
		//}
	//}
	//S += S.transpose();
	//std::cout<<"Symetric Matrix : diagonalizable by orthogonal matrix UU^t=1"<<std::endl;
	//std::cout<<S<<std::endl;
//
	//std::cout<<"The original matrix is preserved (copy inside lapack)"<<std::endl;
	//Lapack<double> Symmetric_copy(S,true,'S');
	//Vector<double> eval;
	//Symmetric_copy.eigensystem(eval,true);
	//Matrix<double> U;
	//Symmetric_copy.set_mat(U);
	//std::cout<<"Eigenvalues"<<std::endl;
	//std::cout<<eval.chop()<<std::endl;
	//std::cout<<"Diagonalization using eigenvector"<<std::endl;
	//std::cout<<((U.transpose()*S*U).diag()).chop()<<std::endl;
	//std::cout<<"Eigenvectors are given in column. Proof :"<<std::endl;
	//Vector<double> v0(S.row());
	//for(unsigned int j(0);j<S.col();j++){
		//for(unsigned int i(0);i<S.row();i++){
			//v0(i) = U(i,j);
		//}
		//std::cout<<((S*v0)-(v0*eval(j))).chop()<<std::endl;
	//}
	//std::cout<<"The original matrix is preserved (explicit copy outside lapack)"<<std::endl;
	//Matrix<double> S_lapack(S);
	//Lapack<double> Symmetric_no_copy(S_lapack,false,'S');
	//eval.set();
	//Symmetric_no_copy.eigensystem(eval,true);
	//Vector<double> vp((S_lapack.transpose()*S*S_lapack).diag().chop());
//
	//std::cout<<"Eigenvalues"<<std::endl;
	//std::cout<<eval<<std::endl;;
	//std::cout<<"Eigenvalues using the passage matrices"<<std::endl;
	//std::cout<<vp<<std::endl;
	///*}*/
	/////*{complex hermitian matrix*/
	////Matrix<std::complex<double> > M(3,3);
	////M(0,0) = std::complex<double> (1,0); 
	////M(1,1) = std::complex<double> (4,0); 
	////M(2,2) = std::complex<double> (5,0); 
	////M(0,1) = std::complex<double> (2,3); 
	////M(1,0) = std::complex<double> (2,-3); 
	////M(0,2) = std::complex<double> (6,4); 
	////M(2,0) = std::complex<double> (6,-4); 
	////M(1,2) = std::complex<double> (-4,-6); 
	////M(2,1) = std::complex<double> (-4,6); 
	////
	////Lapack<std::complex<double> > M_(M,false,'H');
	////Vector<double> EVal2(3);
	////M_.eigensystem(EVal2,true);
	////std::cout<<"eval of a complex matrix"<<std::endl;;
	////std::cout<<EVal2<<std::endl;
	////std::cout<<-7.72113<<" "<<3.4124<<" "<<14.3087<<" (true eigenvalues)" << std::endl;
	/////*}*/
	///*{complex general matrix*/
	//Matrix<std::complex<double> > M(3,3);
	//M(0,0) = std::complex<double> (1,0); 
	//M(1,1) = std::complex<double> (4,0); 
	//M(2,2) = std::complex<double> (5,0); 
	//M(0,1) = std::complex<double> (2,3); 
	//M(1,0) = std::complex<double> (2,3); 
	//M(0,2) = std::complex<double> (10,4); 
	//M(2,0) = std::complex<double> (6,4); 
	//M(1,2) = std::complex<double> (4,1); 
	//M(2,1) = std::complex<double> (-2,6); 
	//
	//Lapack<std::complex<double> > M_(M,true,'G');
	//Matrix<std::complex<double> > EVec;
	//Vector<std::complex<double> > EVal;
	//M_.eigensystem(EVal,&EVec);
	//std::cout<<"eval of a general complex matrix"<<std::endl;
	//std::cout<<EVal<<std::endl;
	//Lapack<std::complex<double> > inv_(EVec,true,'G');
	//inv_.inv();
	//Matrix<std::complex<double> > EVecinv;
	//inv_.set_mat(EVecinv);
//
	//std::cout<<EVal<<std::endl;
	//std::cout<<(EVecinv*M*EVec).chop()<<std::endl;
	///*}*/
	/////*}*/
	///*{lu et det*/
	//Matrix<std::complex<double> > C(3,3);
	//C(0,0)=std::complex<double> (-1.0,3.2);
	//C(0,1)=std::complex<double> (-2.0,1.2);
	//C(0,2)=std::complex<double> (3.0,-7.3);
	//C(1,0)=std::complex<double> (-4.0,1.2);
	//C(1,1)=std::complex<double> (-5.0,1.3);
	//C(1,2)=std::complex<double> (6.0,-6.1);
	//C(2,0)=std::complex<double> (-7.0,7.1);
	//C(2,1)=std::complex<double> (-8.0,8.1);
	//C(2,2)=std::complex<double> (8.0,-1);
	//Lapack<std::complex<double> > C_(&C,false,'G');
	//
	//std::cout<<C<<std::endl;;
	//std::cout<<"det=(160.76,-117.943)="<<C_.det()<<std::endl;
	//
	//Matrix<double> T(3,3);
	//T(0,0)=-1.0;
	//T(0,1)=2.5;
	//T(0,2)=5;
	//T(1,0)=0;
	//T(1,1)=-12;
	//T(1,2)=-145.42;
	//T(2,0)=	7;
	//T(2,1)=54;
	//T(2,2)=47;
	//Lapack<double> Tdet(&T,true,'G');
	//std::cout<<"det=-9413.53="<<Tdet.det()<<std::endl;
	//Lapack<double> Tlu(&T,true,'G');
	//Matrix<double> L(3,3,0.0),U(3,3,0.0);
	//Tlu.lu(L,U);
	//
	//std::cout<<"T"<<std::endl;
	//std::cout<<T<<std::endl;;
	//std::cout<<"L"<<std::endl;
	//std::cout<<L<<std::endl;;
	//std::cout<<"U"<<std::endl;
	//std::cout<<U<<std::endl;;
	//std::cout<<"T=LU"<<std::endl;
	//std::cout<<L*U<<std::endl;;
	///*}*/
	///*{inverse*/
	//Matrix<double> T1(3,3);
	//T1(0,0)=-1.0;
	//T1(0,1)=2.5;
	//T1(0,2)=5;
	//T1(1,0)=0;
	//T1(1,1)=-12;
	//T1(1,2)=-145.42;
	//T1(2,0)=	7;
	//T1(2,1)=54;
	//T1(2,2)=47;
	//Matrix<double> T1inv(T1);
	//
	//std::cout<<"T1 : well defined"<<std::endl;
	//std::cout<<T1<<std::endl;
	//Lapack<double> Over(T1inv,false,'G'); // Tinv va être écrasé
	//double rcond(0.0);	
	//Vector<int> P1(Over.is_singular(rcond));
	//Over.inv(P1);
	//
	//std::cout<<"Tinv"<<std::endl;
	//std::cout<<T1inv<<std::endl;
	//
	//std::cout<<"I"<<std::endl;
	//std::cout<<(T1*T1inv).chop()<<std::endl;;
	//std::cout<<std::endl;
	//
	//Matrix<double> T2(3,3);
	//T2(0,0)=-1.0;
	//T2(0,1)=2.5;
	//T2(0,2)=2;
	//T2(1,0)=0;
	//T2(1,1)=-12;
	//T2(1,2)=-145.42;
	//T2(2,0)=0;
	//T2(2,1)=0;
	//T2(2,2)=1e-200;
	//
	//Lapack<double> Keep(T2,true,'G'); // T2 va être conservé
	//std::cout<<"T2 looks singular"<<std::endl;
	//std::cout<<T2<<std::endl;
	//Vector<int> P2(Keep.is_singular(rcond));
	//Keep.inv(P2);
	//Matrix<double> T2inv_lapack;
	//Keep.set_mat(T2inv_lapack);
	//
	//std::cout<<"get Tinv from lapack --> no inversion has been made"<<std::endl;
	//std::cout<<T2inv_lapack<<std::endl;
	//std::cout<<std::endl;
	//
	//Lapack<double> Keep2(T2,true,'G'); // T2 va être conservé
	//Keep2.inv();
	//std::cout<<"T2 looks singular but will be inverted anyway"<<std::endl;
	//Matrix<double> T2inv_lapack2;
	//Keep2.set_mat(T2inv_lapack2);
//
	//std::cout<<"get T2inv from lapack"<<std::endl;
	//std::cout<<T2inv_lapack2<<std::endl;
	//std::cout<<"result of the inversion of the singular matrix I="<<std::endl;
	//std::cout<<T2inv_lapack2*T2<<std::endl;
	//std::cout<<std::endl;
	//
	//Matrix<std::complex<double> > C(3,3);
	//C(0,0)=std::complex<double>(-1.0,3);
	//C(0,1)=std::complex<double>(2.5,4);
	//C(0,2)=std::complex<double>(2,5);
	//C(1,0)=std::complex<double>(0,0);
	//C(1,1)=std::complex<double>(-12,7);
	//C(1,2)=std::complex<double>(-145.42,7);
	//C(2,0)=std::complex<double>(0,0);
	//C(2,1)=std::complex<double>(0,0);
	//C(2,2)=std::complex<double>(1e-200,2);
	//
	//
	//Lapack<std::complex<double> > inv(C,true,'G');
	//Vector<int> P3(inv.is_singular(rcond));
	//inv.inv(P3);
	//
	//Matrix<std::complex<double> > Cinv;
	//inv.set_mat(Cinv);
	//std::cout<<"complex matrix"<<std::endl;
	//std::cout<<C<<std::endl;
	//std::cout<<"its inverse matrix"<<std::endl;
	//std::cout<<Cinv<<std::endl;
	//std::cout<<"I"<<std::endl;
	//std::cout<<(C*Cinv).chop()<<std::endl<<std::endl;
	////std::cout<<"check that there is no memory leaks"<<std::endl;;
	////for(unsigned int i(0);i<10000000;i++){
	////Lapack<double> Over(&T1,true,'G'); 
	////Over.inv();
	////Matrix<double> Tinv_new(Over.get_mat());
	////std::cout<<(Tinv_new*T1)<<std::endl;;
	////}
	///*}*/
	///*{qr factorisation*/
	//Matrix<double> T(2,5);
	//IOFiles r("QR.dat",false);
	//r>>T;
	//
	//T = T.transpose();
	//Matrix<double>R;
	//Matrix<double>Q;
	//Lapack<double>qr(&T, true, 'G');
	//qr.qr(Q,R,true);
	//Matrix<double>P(qr.get_mat());
	//
	//std::cout<<T.row()<<" "<<T.col()<<std::endl;
	//std::cout<<std::endl;
	//std::cout<<R<<std::endl;
	//std::cout<<std::endl;
	//std::cout<<(T-Q*R*P).chop()<<std::endl;
	///*}*/
	///*{sort*/
	//srand(time(NULL));
	//Vector<double> x(10,1.2);
	//for(unsigned int i(1);i<10;i++){
	//x(i) = rand() % 10 + 1;
	//}
	//Vector<double> x_sorted(x);
	//Vector<unsigned int> index(x_sorted.sort());
	//std::cout<<"original " <<x<<std::endl;
	//std::cout<<"trillée "<<x_sorted<<std::endl;
	//std::cout<<"trillée "<<x.sort(index)<<std::endl;
	///*}*/
	///*{constructors check*/
	//Matrix<double> m1;
	//Matrix<double> m2(2,2) ;
	//Matrix<double> m3(m2);
	//Matrix<double> m4(create_to_check_constructor());
	//m1 = m4;
	//Matrix<double> m5(&m1);
	///*}*/
}

std::complex<double> projection(Matrix<double> const& O, Matrix<std::complex<double> > const& base, unsigned int bra, unsigned int ket){
	unsigned int n(O.row());
	Vector<std::complex<double> > tmp(n,0.0);
	std::complex<double> out(0.0);;
	for(unsigned int i(0);i<n;i++){
		for(unsigned int j(0);j<n;j++){
			tmp(i) += O(i,j)*base(j,ket);
		}
	}
	for(unsigned int i(0);i<n;i++){
		out += tmp(i)*std::conj(base(i,bra));
	}
	return out;
}

Matrix<double> create_to_check_constructor(){
	Matrix<double> m(3,2,2);
	return m;
}
