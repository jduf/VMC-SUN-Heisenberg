#include <armadillo>
#include "State.hpp"
<<<<<<< HEAD
#include "Chrono.hpp"

double energie(unsigned int N_spin, unsigned int N_m, arma::Mat<double> *U);
=======
#include "System.hpp"
#include "Chrono.hpp"

double energie(System *S);
>>>>>>> use-system

int main(){
	Chrono t;
	t.tic();
	for(unsigned int i(0);i<10;i++){
		unsigned int const N_spin(3);
		unsigned int const N_m(6);
<<<<<<< HEAD
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
	}

=======

		System oneD(N_spin, N_m, 1);
		std::cout<<energie(&oneD)<<" ";
	}
>>>>>>> use-system
	t.tac();
	std::cout<<" en "<<t<<std::endl;
}

<<<<<<< HEAD
double energie(unsigned int N_spin, unsigned int N_m, arma::Mat<double> *U){
	unsigned int const N_site(N_spin*N_m); 
	State alpha(N_spin,N_m,U);
	State beta(N_site);
	State tmp(N_site);
=======
double energie(System *S){
	State alpha(S);
	State beta(S->N_site);
	State tmp(S->N_site);
>>>>>>> use-system
	double ratio(0.0), energie(0.0);
	unsigned int i(0),NMC(1e5);
	while(i<NMC){
		//std::cout<<"tmp ";
		tmp = alpha.swap();
		ratio = (tmp.Det() * tmp.Det()) / (alpha.Det() * alpha.Det());
		if(ratio>1 || rand()/RAND_MAX <ratio){
<<<<<<< HEAD
			alpha = tmp;
			i++;
			for(unsigned int j(0); j<N_site;j++){
				beta = alpha.swap(j,(j+1) % N_site);
=======
			i++;
			//std::cout<<"alpha ";
			alpha = tmp;
			for(unsigned int j(0);j<S->N_site;j++){
				//std::cout<<"beta ";
				beta = alpha.swap(j,S->nts[j]);
>>>>>>> use-system
				energie += beta.Det()/alpha.Det();
			}
		} 
		//std::cout<<std::endl;
	}
	// sign(permutation) => ratio det tjrs - ??? 
	return -energie/(S->N_site * NMC);
}


