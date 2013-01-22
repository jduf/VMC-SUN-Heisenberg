#include "State.hpp"
#include "System.hpp"
#include "Chrono.hpp"

double energie(System *S);

int main(){
	Chrono t;
	t.tic();
	unsigned int const N_spin(3);
	unsigned int const N_m(2);
	System oneD(N_spin, N_m, 1);	
	std::cout<<energie(&oneD)<<" ";
	t.tac();
	std::cout<<"en "<<t<<" seconde(s)"<<std::endl;
}

double energie(System *S){
	State alpha(S);
	State beta(S->N_site,S->N_spin);
	State tmp(S->N_site,S->N_spin);

	double ratio(0.0), energie(0.0);
	unsigned int i(0),NMC(30);
	while(i<NMC){
		tmp = alpha.swap();
		ratio = (tmp.Det() * tmp.Det()) / (alpha.Det() * alpha.Det());
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			i++;
			alpha = tmp;
			for(unsigned int j(0);j<S->N_site;j++){
				for(unsigned int d(0);d<S->dim;d++){
					beta = alpha.swap(j,S->nts[S->dim*j+d]);
					energie += beta.Det()/alpha.Det();
				}
			}
		} else { std::cout<<"pas "<<i<<std::endl; }
	}
	return -energie/(S->N_site * NMC); // sign(permutation) => ratio det tjrs - ??? 
}


