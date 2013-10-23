/**
  @file test.cpp 
  */

#include "Gnuplot.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"
#include "Lapack.hpp"
#include "Read.hpp"
#include "Write.hpp"

#include <complex>
#include <stdlib.h>
#include <time.h>

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

int main(){
	/*operateurs*/
	///*{*/
	//unsigned int N1_row(3);
	//unsigned int N1_col(2);
	//unsigned int N2_row(2);
	//unsigned int N2_col(6);
	//Matrix<double> m1(N1_row,N1_col);
	//Matrix<double> m2(N2_row,N2_col,2);
	//
	//m1(0,0) = 1;
	//m1(0,1) = 2;
	//m1(1,0) = 3;
	//m1(1,1) = 4;
	//m1(2,0) = 5;
	//m1(2,1) = 6;
	//Matrix<double> m3(m1);
	//
	//std::cout<<"m1"<<std::endl;
	//std::cout<<m1<<std::endl;
	//std::cout<<"m2"<<std::endl;
	//std::cout<<m2<<std::endl;
	//std::cout<<"m3"<<std::endl;
	//std::cout<<m3<<std::endl;
	//m1 = m1*m2;
	//std::cout<<"m1=m1*m2"<<std::endl;
	//std::cout<<m1<<std::endl;
	//m2=m2.transpose();
	//std::cout<<"m2=m2.tanspose"<<std::endl;
	//std::cout<<m2<<std::endl;
	//std::cout<<"m1*m2"<<std::endl;
	//std::cout<<m1*m2<<std::endl;
	//std::cout<<"m3-(m1-m2)"<<std::endl;
	//std::cout<<m3-(m1*m2)<<std::endl;
	//
	//Matrix<double> m4(1,3,3.1);
	//Matrix<double> m5(3,1,11);
	//
	//std::cout<<"m4*m5"<<std::endl;
	//std::cout<<m5*m4<<std::endl;
	//m4 *= m1;
	//std::cout<<m4<<std::endl;
	//m1.print_mathematica();
	//std::cout<<m1<<std::endl;
	// /*}*/
	/*eigenvalue*/
	/*{*/
	/*Matrix initialization*/
	/*{*/

	//unsigned int N_site(11);
	//Matrix<double> H(N_site,N_site,0.0);
	//H(0,1)=-1.0;
	//H(0,N_site-1)=1.0;
	//for(unsigned int i(1); i< N_site-1; i++){
		//H(i,i-1) = -1.0;
		//H(i,i+1) = -1.0;
	//}
	//H(N_site-1,0)=1.0;
	//H(N_site-1,N_site-2)=-1.0;
	/*}*/

	/*projection <ket|O|bra>*/
	///*{*/
	//Matrix<double> P(N_site,N_site,0.0);
	//P(N_site-1,0)=-1.0;
	//for(unsigned int i(0); i< N_site-1; i++){
	//P(i,i+1) = 1.0;
	//}
		//
	//Matrix<double> HP(H+P);
	//Lapack<double> HP_(&HP,true,'G');
	//Matrix<std::complex<double> > eve;
	//Vector<std::complex<double> > eva;
	//HP_.eigensystem(&eva,&eve);
	//std::cout<<eva<<std::endl;
	//
	//Vector<double> k(N_site);
	//Vector<double> E(N_site);
	//for(unsigned int i(0);i<N_site;i++){
		//k(i) = log(projection(P,eve,i,i)).imag();
		//E(i) = projection(H,eve,i,i).real();
	//}
	//
	//Gnuplot gp("spectrum","1D");
	//gp.save_data("spectrum",k,E);
	//gp.save_code();
	///*}*/
	///*{*/
	//Lapack<double> U_(&H,true,'S');
	//Matrix<double> EVal;
	//U_.eigensystem(&EVal,true);
	//Matrix<double> U(U_.get_mat());
	//
	//std::cout<<((U.transpose()*H*U).diag()).transpose().chop()<<std::endl;
	//std::cout<<EVal.transpose().chop()<<std::endl;
	//Matrix<double> v0(N_site,1);
	//for(unsigned int i(0);i<N_site;i++){
	//v0(i) = U(i,0);
	//}
	//std::cout<<((H*v0)-(v0*EVal(0))).transpose().chop()<<std::endl;
	///*}*/
	/*check how eigenvectors are stored*/
	///*{*/
	//std::cout<<"Original matrix"<<std::endl;
	//std::cout<<H<<std::endl;;
	//Matrix<double> T(H);
	//Vector<double> eval;
	//Lapack<double> ES(&T,false,'S');
	//ES.eigensystem(&eval,true);
//
	//Matrix<double> Tinv(T); // T is now the passage matrix
	//Lapack<double> Tinv_(&Tinv,false,'G'); 
	//Tinv_.inv();
	//Matrix<double> vp(Tinv*H*T);
//
	//std::cout<<"Eigenvalues"<<std::endl;
	//std::cout<<eval.chop()<<std::endl;;
	//std::cout<<"Eigenvalues using the passage matrices"<<std::endl;
	//std::cout<<vp.diag().chop()<<std::endl;
	///*}*/
	///*{*/
	//std::cout<<"Original matrix"<<std::endl;
	//std::cout<<H<<std::endl;;
	//Matrix<double> T(H);
	//Matrix<double> EValOW(N_site,1);
	//Matrix<double> EValNOW(N_site,1);
//
	//Lapack<double> EVecNOW(&T,true,'S');
	//Lapack<double> EVecOW(&T,false,'S');
	//EVecNOW.eigensystem(&EValNOW,true);
	//EVecOW.eigensystem(&EValOW,true);
	//EValNOW.chop();
	//EValOW.chop();
//
	//Matrix<double> Tinv(T); // T is now the passage matrix
	//Lapack<double> Tinv_(&Tinv,false,'G'); 
	//Tinv_.inv();
	//Matrix<double> vp(Tinv*H*T);
//
	//std::cout<<"Eigenvalues without overwriting the original matrix"<<std::endl;
	//std::cout<<EValNOW.chop()<<std::endl;;
	//std::cout<<"Eigenvalues with overwriting the original matrix"<<std::endl;
	//std::cout<<EValOW.chop()<<std::endl;;
	//std::cout<<"Eigenvalues using the passage matrices"<<std::endl;
	//std::cout<<vp.diag().chop()<<std::endl;
//
	//T.chop();
	//std::cout<<"Passage matrix : Tinv.H.T = eigenvalues"<<std::endl;
	//std::cout<<T<<std::endl;;
	//std::cout<<std::endl;;
	///*}*/
	/*complex hermitian matrix*/
	///*{*/
	//Matrix<std::complex<double> > M(3,3);
	//M(0,0) = std::complex<double> (1,0); 
	//M(1,1) = std::complex<double> (4,0); 
	//M(2,2) = std::complex<double> (5,0); 
	//M(0,1) = std::complex<double> (2,3); 
	//M(1,0) = std::complex<double> (2,-3); 
	//M(0,2) = std::complex<double> (6,4); 
	//M(2,0) = std::complex<double> (6,-4); 
	//M(1,2) = std::complex<double> (-4,-6); 
	//M(2,1) = std::complex<double> (-4,6); 
//
	//Lapack<std::complex<double> > M_(&M,false,'H');
	//Vector<double> EVal2(3);
	//M_.eigensystem(&EVal2,true);
	//std::cout<<"eval of a complex matrix"<<std::endl;;
	//std::cout<<EVal2<<std::endl;
	//std::cout<<-7.72113<<" "<<3.4124<<" "<<14.3087<<" (true eigenvalues)" << std::endl;
	///*}*/
	/*complex general matrix*/
	/*{*/
	Matrix<std::complex<double> > M(3,3);
	M(0,0) = std::complex<double> (1,0); 
	M(1,1) = std::complex<double> (4,0); 
	M(2,2) = std::complex<double> (5,0); 
	M(0,1) = std::complex<double> (2,3); 
	M(1,0) = std::complex<double> (2,3); 
	M(0,2) = std::complex<double> (10,4); 
	M(2,0) = std::complex<double> (6,4); 
	M(1,2) = std::complex<double> (4,1); 
	M(2,1) = std::complex<double> (-2,6); 

	Lapack<std::complex<double> > M_(&M,false,'G');
	Matrix<std::complex<double> > EVec(3,3);
	Vector<std::complex<double> > EVal(3);
	M_.eigensystem(&EVal,&EVec);
	std::cout<<"eval of a genral complex matrix"<<std::endl;;
	std::cout<<EVal<<std::endl;
	/*}*/
	/*}*/
	/*lu et det*/
	///*{*/
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
	//inverse
	///*{*/
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
	//Lapack<double> Over(&T1inv,false,'G'); // Tinv va être écrasé
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
	//Lapack<double> Keep(&T2,true,'G'); // T2 va être conservé
	//std::cout<<"T2 looks singular"<<std::endl;
	//std::cout<<T2<<std::endl;
	//Vector<int> P2(Keep.is_singular(rcond));
	//Keep.inv(P2);
	//Matrix<double> T2inv_lapack(Keep.get_mat());
//
	//std::cout<<"get Tinv from lapack --> no inversion has been made"<<std::endl;
	//std::cout<<T2inv_lapack<<std::endl;
	//std::cout<<std::endl;
//
	//Lapack<double> Keep2(&T2,true,'G'); // T2 va être conservé
	//Keep2.inv();
	//std::cout<<"T2 looks singular but will be inverted anyway"<<std::endl;
	//Matrix<double> T2inv_lapack2(Keep2.get_mat());
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
	//Lapack<std::complex<double> > inv(&C,true,'G');
	//Vector<int> P3(inv.is_singular(rcond));
	//inv.inv(P3);
//
	//Matrix<std::complex<double> > Cinv(inv.get_mat());
	//std::cout<<"complex matrix"<<std::endl;
	//std::cout<<C<<std::endl;
	//std::cout<<"its inverse matrix"<<std::endl;
	//std::cout<<Cinv<<std::endl;
	//std::cout<<"I"<<std::endl;
	//std::cout<<(C*Cinv).chop()<<std::endl;
//
	////std::cout<<"check that there is no memory leaks"<<std::endl;;
	////for(unsigned int i(0);i<10000000;i++){
		////Lapack<double> Over(&T1,true,'G'); 
		////Over.inv();
		////Matrix<double> Tinv_new(Over.get_mat());
		////std::cout<<(Tinv_new*T1)<<std::endl;;
	////}
	///*}*/
	//qr factorisation
	///*{*/
	//Matrix<double> T;
	//Read r("../../SUN/dev/src/sim/test_something.jdbin");
	//r>>T;
	//
	////std::cout<<T<<std::endl;
	////std::cout<<T*P<<std::endl;
	////Matrix<double> A(T.row(),T.row());
	////Matrix<double> tmp(T*P);
	////for(unsigned int i(0);i<A.row();i++){
	////for(unsigned int j(0);j<A.row();j++){
	////A(i,j) = tmp(i,j);
	////}
	////}
	//
	////Lapack<double> tmp_(&A,true,'G');
	////std::cout<<tmp_.det()<<std::endl;
	//
	//
	//T = T.transpose();
	//Matrix<double>R;
	//Matrix<double>Q;
	//Lapack<double>qr(&T, true, 'G');
	//qr.qr(Q,R,true);
	//Matrix<double>P(qr.get_mat());
	//T.chop();
	//Q.chop();
	//R.chop();
	//T.chop();
	//
	////std::cout<<Q<<std::endl;
	////std::cout<<R<<std::endl;
	//std::cout<<P<<std::endl;
	//std::cout<<T<<std::endl;
	//std::cout<<Q*R*P<<std::endl;
	////Matrix<double> ouf(T.col(),T.col());
	////Matrix<double> QR(Q*R);
	////for(unsigned int i(0);i<ouf.row();i++){
	////for(unsigned int j(0);j<ouf.col();j++){
	////ouf(i,j) = QR(i,j);
	////}
	////}
	////Lapack<double> det_(&ouf,false,'G');
	////std::cout<<det_.det()<<std::endl;
	//
	//Write wt("T.dat");
	//wt<<T;
	//Write wqr("QR.dat");
	//wqr<<Q*R*P;
	//Write wq("Q.dat");
	//wq<<Q;
	//Write wr("R.dat");
	//wr<<R;
	///*}*/
	/*sort*/
	///*{*/
	//srand(time(NULL));
	//Vector<double> x(10,1.2);
	//for(unsigned int i(1);i<10;i++){
		//x(i) = rand() % 10 + 1;
	//}
	//Vector<double> x_sorted(x);
	//Vector<unsigned int> index(x_sorted.sort());
	//std::cout<<x<<std::endl;
	//std::cout<<std::endl;
	//std::cout<<x_sorted<<std::endl;
	//std::cout<<std::endl;
	//std::cout<<x.sort(index)<<std::endl;
//
	//for(unsigned int i(0);i<10;i++){
	//std::cout<<x[i] <<" ";
	//}
	//std::cout<<x.sort().transpose()<<std::endl;
	///*}*/
	
	/*constructors check*/
	///*{*/
	//Matrix<double> m1;
	//Matrix<double> m2(2,2) ;
	//Matrix<double> m3(m2);
	//Matrix<double> m4(create_to_check_constructor());
	//m1 = m4;
	//Matrix<double> m5(&m1);
	///*}*/
}
