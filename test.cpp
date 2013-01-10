#include <armadillo>
#include "State.hpp"
#include "System.hpp"

double energie(System *S);

int main(){
	unsigned int const N_spin(3);
	unsigned int const N_m(5);

	System oneD(N_spin, N_m, 1);

	std::cout<<energie(&oneD)<<std::endl;
}

double energie(System *S){
	State alpha(S);
	State beta;
	State tmp;
	double ratio(0.0), energie(0.0);
	unsigned int i(0),NMC(1e6);
	while(i<NMC){
		tmp = alpha.swap();
		ratio = (tmp.Det() * tmp.Det()) / (alpha.Det() * alpha.Det());
		if(ratio>1 || rand()/RAND_MAX <ratio){
			i++;
			alpha = tmp;
			for(unsigned int j(0);j<S->N_site;j++){
				beta = alpha.swap(j,S->nts[j]);
				energie += beta.Det()/alpha.Det();
			}
		} 
	}
	// sign(permutation) => ratio det tjrs - ??? 
	return -energie/(S->N_site * NMC);
}


