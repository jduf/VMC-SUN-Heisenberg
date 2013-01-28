#include "Save.hpp"
#include "State.hpp"
#include "System.hpp"
#include "Chrono.hpp"

double energie(System *S);

int main(){
	//Chrono t;
	//t.tic();
	unsigned int const N_spin(3);
	unsigned int const N_m(4);
	System oneD(N_spin, N_m, 1);	
	//std::cout<<energie(&oneD)<<" ";
	//t.tac();
	//std::cout<<"en "<<t<<" seconde(s)"<<std::endl;
	
	std::cout<<"create alpha"<<std::endl;
	State alpha(&oneD);
	alpha.print();
	State tmp(alpha.swap());
	tmp.print();
	std::cout<<tmp/alpha<<std::endl;

	//t.tic();
	//unsigned int const N_spin(4);
	//unsigned int const N_m(4);
	//System two(N_spin, N_m, 2);	
	//std::cout<<energie(&two)<<" ";
	//t.tac();
	//std::cout<<"en "<<t<<" seconde(s)"<<std::endl;
}

double energie(System *S){
	State alpha(S);
	double ratio(0.0), energie(0.0);
	unsigned int i(0),NMC(1e4);
	while(i<NMC){
		State tmp(alpha.swap());
		ratio = tmp/alpha;
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			i++;
			alpha = tmp;
			for(unsigned int j(0);j<S->N_site;j++){
				for(unsigned int d(0);d<S->dim;d++){
					State beta(alpha.swap(j,S->nts[S->dim*j+d]));
					energie += beta/alpha;
				}
			}
		}
	}
	return -energie/(S->N_site * NMC); // sign(permutation) => ratio det tjrs - ??? 
}


