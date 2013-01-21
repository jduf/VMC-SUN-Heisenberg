#include "State.hpp"
#include "System.hpp"
#include "Chrono.hpp"

double energie(System *S);

int main(){
	//Chrono t;
	//t.tic();
	unsigned int const N_spin(3);
	unsigned int const N_m(2);
	System oneD(N_spin, N_m, 1);
	std::cout<<energie(&oneD)<<" ";
	//t.tac();
	//std::cout<<"en "<<t<<std::endl;
}

double energie(System *S){
	std::cout<<"initialisation states S"<<std::endl;
	State alpha(S);
	std::cout<<std::endl;
	std::cout<<"initialisation states beta"<<std::endl;
	State beta(S->N_site,S->N_spin);
	std::cout<<std::endl;
	std::cout<<"initialisation states tmp"<<std::endl;
	State tmp(S->N_site,S->N_spin);
	std::cout<<std::endl;

	double ratio(0.0), energie(0.0);
	unsigned int i(0),NMC(6);
	while(i<NMC){
		std::cout<<"swap alpha"<<std::endl;
		tmp = alpha.swap();
		std::cout<<"det computed"<<std::endl;
		ratio = (tmp.Det() * tmp.Det()) / (alpha.Det() * alpha.Det());
		if(ratio>1 || rand()/RAND_MAX <ratio){
			i++;
			std::cout<<"garde alpha"<<std::endl;
			alpha = tmp;
			for(unsigned int j(0);j<S->N_site;j++){
				for(unsigned int d(0);d<S->dim;d++){
					beta = alpha.swap(j,S->nts[S->dim*j+d]);
					energie += beta.Det()/alpha.Det();
				}
			}
		} 
	}
	return -energie/(S->N_site * NMC); // sign(permutation) => ratio det tjrs - ??? 
}


