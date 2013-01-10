#include "System.hpp"

System::System(unsigned int N_spin, unsigned int N_m, unsigned int dim):
	N_spin(N_spin),
	N_m(N_m),
	N_site(N_m*N_spin),
	nts(new unsigned int[N_spin*N_m*dim]),
	U()
{
	create_U(dim);
	create_nts(dim);
}

System::~System(){
	delete nts;
}

void System::create_U(unsigned int dim){
	arma::Mat<double> T(arma::zeros(N_site,N_site));
	if(dim==1){
		T(0,1)=-1.0;
		T(0,N_site-1)=1.0;
		for(unsigned int i(1); i< N_site-1; i++){
			T(i,i-1) = -1.0;
			T(i,i+1) = -1.0;
		}
		T(N_site-1,0)=1.0;
		T(N_site-1,N_site-2)=-1.0;

	}
	
	arma::Col<double> EVal;
	arma::eig_sym(EVal,U,T);
	U.print();
}

void System::create_nts(unsigned int dim){
	if(dim==1){
		for(unsigned int i(0);i<N_site;i++){
			nts[i] = (i+1) % N_site;
		}
	}
}
