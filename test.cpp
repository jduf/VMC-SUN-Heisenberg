#include "Save.hpp"
#include "System.hpp"
#include "Chrono.hpp"

double energie(System& S);

int main(){
	unsigned int N_spin(3),N_m(12),N_MC(1e4);
	//std::cout<<"N_spin="<<std::flush;
	//get_N(N_spin);
	//std::cout<<"N_m="<<std::flush;
	//get_N(N_m);
	//std::cout<<"N_MC="<<std::flush;
	//get_N(N_MC);
	//std::cout<<N_MC<<std::endl;

	System S(N_spin, N_m, 2);

	Chrono t;
	t.tic();
	std::cout<<energie(S)<<std::endl;
	t.tac();
	std::cerr<<t<<" seconde(s)"<<std::endl;
}

double energie(System& S){
	double ratio(0.0), energie(0.0);
	unsigned int i(0);
	Save steps("analyse-2d.dat");
	while(i<N_MC){
		S.swap();
		ratio = S.compute_ratio();
		ratio *= ratio;
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			i++;
			S.update_state();
			steps << S << " " << S.det() <<Save::endl;
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


