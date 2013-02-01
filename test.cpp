#include "Save.hpp"
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
	std::cout<<energie(&oneD)<<std::endl;
	t.tac();
	std::cerr<<t<<" seconde(s)"<<std::endl;
}

double energie(System *S){
	State alpha(S);

	double ratio(0.0), energie(0.0);
	unsigned int i(0),NMC(40);
	while(i<NMC){
		State tmp(alpha.swap());
		ratio = (tmp.Det() * tmp.Det()) / (alpha.Det() * alpha.Det());
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			i++;
			alpha = tmp;
			for(unsigned int j(0);j<S->N_site;j++){
				for(unsigned int d(0);d<S->dim;d++){
					State beta(alpha.swap(j,S->nts[S->dim*j+d]));
					energie += beta.Det()/alpha.Det();
				}
			}
		}
	}
	return -energie/(S->N_site * NMC); // sign(permutation) => ratio det tjrs - ??? 
}


