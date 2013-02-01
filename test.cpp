#include "Save.hpp"
#include "System.hpp"
#include "Chrono.hpp"

double energie(System& S);

int main(){
	unsigned int const N_spin(3);
	unsigned int const N_m(2);
	System oneD(N_spin, N_m, 1);	

	Chrono t;
	t.tic();
	std::cout<<energie(oneD)<<std::endl;
	t.tac();
	std::cerr<<t<<" seconde(s)"<<std::endl;
}

double energie(System& S){
	double ratio(0.0), energie(0.0);
	unsigned int i(0),NMC(40);
	while(i<NMC){
		S.swap();
		ratio = S.compute_ratio();
		ratio *= ratio;
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			i++;
			S.update_state();
			//S.print();
			for(unsigned int j(0);j<S.N_site;j++){
				for(unsigned int d(0);d<S.dim;d++){
					S.swap(j,S.nts[S.dim*j+d]);
					energie += S.compute_ratio();
				}
			}
		}
	}
	return -energie/(S.N_site * NMC); // sign(permutation) => ratio det tjrs - ??? 
}


