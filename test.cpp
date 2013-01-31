#include "Save.hpp"
#include "State.hpp"
#include "System.hpp"
#include "Chrono.hpp"

double energie(System *S);

int main(){
	unsigned int const N_spin(3);
	unsigned int const N_m(4);
	System oneD(N_spin, N_m, 1);	

	//Chrono t;
	//t.tic();
	//std::cout<<energie(&oneD)<<std::endl;
	//t.tac();
	//std::cerr<<t<<" seconde(s)"<<std::endl;
	
	State alpha(&oneD,true);
	State tmp(&oneD,false);

	//alpha.print();
	tmp = alpha.swap();
	std::cout<<tmp/alpha<<std::endl;
	alpha=tmp;
	alpha.print();
}

double energie(System *S){
	State alpha(S,true);
	State tmp(S,false);
	double ratio(0.0), energie(0.0);
	unsigned int i(0),NMC(40);
	while(i<NMC){
		tmp = alpha.swap();
		ratio = tmp/alpha;
		ratio *= ratio;
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			i++;
			alpha = tmp;
			for(unsigned int j(0);j<S->N_site;j++){
				for(unsigned int d(0);d<S->dim;d++){
					tmp = alpha.swap(j,S->nts[S->dim*j+d]);
					energie += tmp/alpha;
				}
			}
		}
	}
	return -energie/(S->N_site * NMC); // sign(permutation) => ratio det tjrs - ??? 
}


