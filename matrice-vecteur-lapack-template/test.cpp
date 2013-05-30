/**
@file test.cpp 
*/

#include "Vecteur.hpp"
#include "Matrice.hpp"
#include "Lapack.hpp"

#include <complex>

int main(){
	/*operateurs*/
	///*{*/
	//unsigned int N(2);
	//Matrice<double> m1(N);
	//Matrice<double> m2(N,2);
	////
	//m1(0,0) = 1;
	//m1(0,1) = 2;
	//m1(1,0) = 3;
	//m1(1,1) = 4;
	//m2(0,1) = 5;
	//m2(1,0) = 6;
	//m2(1,1) = 7;
//
	//std::cout<<"m1"<<std::endl;
	//std::cout<<m1<<std::endl;
	//std::cout<<"m2"<<std::endl;
	//std::cout<<m2<<std::endl;
	//std::cout<<"m1*m2"<<std::endl;
	//std::cout<<m1*m2<<std::endl;
	//m1 *= m2;
	//std::cout<<m1<<std::endl;
	//
	//std::cout<<"m2-(m1*m2)"<<std::endl;
	//std::cout<<m2-m1<<std::endl;
	//m2-=m1;
	//std::cout<<"(m2-(m1*m2))-(m1*m2)"<<std::endl;
	//std::cout<<m2-m1<<std::endl;
//
	//Vecteur<double> v1(3,3);
	//std::cout<<"v1"<<std::endl;
	//std::cout<<v1<<std::endl;
	//Vecteur<double> v2(3);
	//v2(0) = 0.0;
	//v2(1) = v1(1);
	//v2(2) = 5;
	//std::cout<<"v2"<<std::endl;
	//std::cout<<v2<<std::endl;
//
	//std::cout<<"v1^v2"<<std::endl;;
	//std::cout<<(v1^v2)<<std::endl;;
	//std::cout<<"v1*v2"<<std::endl;
	//std::cout<<v1*v2<<std::endl;
//
	//std::cout<<"v1*2"<<std::endl;
	//std::cout<<v1*2<<std::endl;
	//v1 *= 2;
	//std::cout<<v1<<std::endl;
	//std::cout<<"2*v2"<<std::endl;
	//std::cout<<2.0*v2<<std::endl;
	//
	//std::cout<<"v2/3.0"<<std::endl;
	//std::cout<<v2/3.0<<std::endl;
	//v2/=3.0;
	//std::cout<<v2<<std::endl;
//
	//Vecteur<double> vtest(0);
	//vtest = Vecteur<double>(3,1);
	//std::cout<<"test"<<vtest<<std::endl;
//
	//std::complex<double> c1(1,2);
	//std::complex<double> c2(4,5);
	//std::complex<double> c3(5,1);
	//Vecteur<std::complex<double> > v3(3,c1);
	//Vecteur<std::complex<double> > v4(3);
//
	//v4(0) = v3(0);
	//v4(1) = c3;
	//v4(2) = 2;
	//
	//std::cout<<"v3"<<std::endl;
	//std::cout<<v3<<std::endl;
	//std::cout<<"v4"<<std::endl;
	//std::cout<<v4<<std::endl;
	//std::cout<<"v3*v4"<<std::endl;
	//std::cout<<v3*v4<<std::endl;
	//std::cout<<"v3*(4,5)"<<std::endl;
	//std::cout<<v3*c2<<std::endl;
	//std::cout<<"(4,5)*v3"<<std::endl;
	//std::cout<<c2*v3<<std::endl;
//
	// /*}*/
	/*eigenvalue*/
	/*{*/
	unsigned int N_site(6);
	Matrice<double> H(N_site,0.0);
	H(0,1)=-1.0;
	H(0,N_site-1)=1.0;
	for(unsigned int i(1); i< N_site-1; i++){
		H(i,i-1) = -1.0;
		H(i,i+1) = -1.0;
	}
	H(N_site-1,0)=1.0;
	H(N_site-1,N_site-2)=-1.0;
	std::cout<<H<<std::endl;;
	Matrice<double> T(H);
	Vecteur<double> EVal(N_site);
	Lapack<double> EVecTest(T,'S');
	EVecTest.eigensystem(EVal);
	EVal.chop();
	std::cout<<EVal<<std::endl;;
	Lapack<double> EVec(T.ptr(),T.size(),'S');
	EVec.eigensystem(EVal,true);
	EVal.chop();
	std::cout<<EVal<<std::endl;;
	T.chop();
	std::cout<<T<<std::endl;;
	
	Matrice<double> Tinv(T);
	Lapack<double> Tinv_(Tinv.ptr(),Tinv.size(),'G');
	Tinv_.inv();
	
	Matrice<double> vp(Tinv*H*T);
	std::cout<<vp.diag()<<std::endl;

	//Matrice<std::complex<double> > M(3,0);
	//M(0,0) = std::complex<double> (1,0); 
	//M(0,1) = std::complex<double> (2,3); 
	//M(0,2) = std::complex<double> (6,4); 
	//M(1,0) = std::complex<double> (2,-3); 
	//M(1,1) = std::complex<double> (4,0); 
	//M(1,2) = std::complex<double> (-4,-6); 
	//M(2,0) = std::complex<double> (6,-4); 
	//M(2,1) = std::complex<double> (-4,6); 
	//M(2,2) = std::complex<double> (5,0); 
//
	//Lapack<std::complex<double> > M_(M.ptr(),M.size(),'H');
	//Vecteur<double> EVal2(3);
	//M_.eigensystem(EVal2,true);
	//std::cout<<EVal2<<std::endl;;
	//std::cout<<"true eigenvalue" << -7.72113<<" "<<3.4124<<" "<<14.3087<<std::endl;
	//std::cout<<M<<std::endl;;
	//std::cout<<M.trans_conj()<<std::endl;;
	/*}*/
	//lu et det
	///*{*/
	//unsigned int N_site(3);
	//Matrice<double> T(N_site);
	//T(0,0)=-1.0;
	//T(0,1)=2.5;
	//T(0,2)=5;
	//T(1,0)=0;
	//T(1,1)=-12;
	//T(1,2)=-145.42;
	//T(2,0)=	7;
	//T(2,1)=54;
	//T(2,2)=47;
	//T<<std::endl;;
	//std::cout<<"det="<<-9413.53<<std::endl;	
	//Lapack<double> T_(T,'G');
	//
	//std::cout<<T_.det()<<std::endl;
//
	//Matrice<std::complex<double> > C(N_site);
	//C(0,0)=std::complex<double> (-1.0,3.2);
	//C(0,1)=std::complex<double> (-2.0,1.2);
	//C(0,2)=std::complex<double> (3.0,-7.3);
	//C(1,0)=std::complex<double> (-4.0,1.2);
	//C(1,1)=std::complex<double> (-5.0,1.3);
	//C(1,2)=std::complex<double> (6.0,-6.1);
	//C(2,0)=std::complex<double> (-7.0,7.1);
	//C(2,1)=std::complex<double> (-8.0,8.1);
	//C(2,2)=std::complex<double> (8.0,-1);
	//C<<std::endl;;
	//std::cout<<"det=160.76-117.943 i"<<std::endl;	
	//Lapack<std::complex<double> > C_(C,'G');
	//
	//std::cout<<C_.det()<<std::endl;
//
	////Matrice L(N_site,0.0),U(N_site,0.0);
	////T_.lu(L,U);
	////std::cout<<"L"<<std::endl;
	////L<<std::endl;;
	////std::cout<<"U"<<std::endl;
	////U<<std::endl;;
	////std::cout<<"LU"<<std::endl;
	////(L*U)<<std::endl;;
	////
	////Matrice A(2);
	////A(0,0) = 1.2;
	////A(0,1) = 1.0;
	////A(1,0) = 1.4;
	////A(1,1) = 2;
////
	////Lapack A_(A,'G');
	////std::cout<<A_.det()<<std::endl;
//
//
	///*}*/
	//inverse
	///*{*/
	//unsigned int N_site(3);
	//Matrice T(N_site);
	//T(0,0)=-1.0;
	//T(0,1)=2.5;
	//T(0,2)=5;
	//T(1,0)=0;
	//T(1,1)=-12;
	//T(1,2)=-145.42;
	//T(2,0)=	7;
	//T(2,1)=54;
	//T(2,2)=47;
	//T<<std::endl;;
	//std::cout<<"det="<<-9413.53<<std::endl;	
//
	//Lapack Keep(T,'G');
	//Keep.inv();
//
	//Matrice Tinv(Keep.LapToMat());
	//Tinv<<std::endl;;
//
	//(T*Tinv)<<std::endl;;
	//(Tinv*T)<<std::endl;;
//
	//Lapack Over(T.ptr(),T.size(),'G'); // T va être écrasé
	//
	//T<<std::endl;;
//
	///*}*/
}



