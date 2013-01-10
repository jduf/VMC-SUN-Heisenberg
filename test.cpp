#include <armadillo>
#include "State.hpp"
#include "System.hpp"
#include "Chrono.hpp"

double energie(System *S);

int main(){
	Chrono t;
	t.tic();
	unsigned int const N_spin(3);
	unsigned int const N_m(6);

	System oneD(N_spin, N_m, 1);
	for(unsigned int i(0);i<10;i++){
		std::cout<<energie(&oneD)<<" ";
	}
	t.tac();
	std::cout<<" en "<<t<<std::endl;
}

double energie(System *S){
	State alpha(S);
	State beta(S->N_site);
	State tmp(S->N_site);
	double ratio(0.0), energie(0.0);
	unsigned int i(0),NMC(1e5);
	while(i<NMC){
		//std::cout<<"tmp ";
		tmp = alpha.swap();
		ratio = (tmp.Det() * tmp.Det()) / (alpha.Det() * alpha.Det());
		if(ratio>1 || rand()/RAND_MAX <ratio){
			i++;
			//std::cout<<"alpha ";
			alpha = tmp;
			for(unsigned int j(0);j<S->N_site;j++){
				//std::cout<<"beta ";
				beta = alpha.swap(j,S->nts[j]);
				energie += beta.Det()/alpha.Det();
			}
		} 
		//std::cout<<std::endl;
	}
	// sign(permutation) => ratio det tjrs - ??? 
	return -energie/(S->N_site * NMC);
}


