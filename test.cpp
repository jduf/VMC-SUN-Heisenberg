#include <armadillo>
#include "State.hpp"
#include "Chrono.hpp"

double energie(unsigned int N_spin, unsigned int N_m, arma::Mat<double> *U);

int main(){
	Chrono t;
	t.tic();
	unsigned int const N_spin(3);
	unsigned int const N_m(6);
	unsigned int const N_site(N_spin*N_m);

	arma::Mat<double> U(arma::zeros(N_site,N_site));
	U(0,1)=-1.0;
	U(0,N_site-1)=1.0;
	for(unsigned int i(1); i< N_site-1; i++){
		U(i,i-1) = -1.0;
		U(i,i+1) = -1.0;
	}
	U(N_site-1,0)=1.0;
	U(N_site-1,N_site-2)=-1.0;
	
	arma::Col<double> EVal;
	arma::Mat<double> EVec;
	arma::eig_sym(EVal,EVec,U);
	
	std::cout<<energie(N_spin,N_m,&EVec);
	t.tac();
	std::cout<<" en "<<t<<std::endl;
}

double energie(unsigned int N_spin, unsigned int N_m, arma::Mat<double> *U){
	unsigned int const N_site(N_spin*N_m); 
	State alpha(N_spin,N_m,U);
	State beta(N_site);
	State tmp(N_site);
	double ratio(0.0), energie(0.0);
	unsigned int i(0),NMC(1e3);
	while(i<NMC){
		tmp = alpha.swap();
		ratio = (tmp.Det() * tmp.Det()) / (alpha.Det() * alpha.Det());
		if(ratio>1 || rand()/RAND_MAX <ratio){
			alpha = tmp;
			i++;
			for(unsigned int j(0); j<N_site;j++){
				beta = alpha.swap(j,(j+1) % N_site);
				energie += beta.Det()/alpha.Det();
			}
		} 
	}
	// sign(permutation) => ratio det tjrs - ??? 
	return -energie/(N_m * N_spin * NMC);
}


