#include "Save.hpp"
#include "Chrono.hpp"
#include "Parseur.hpp"
#include "System.hpp"

#include<string>
#include<sstream>
#include<iostream>

double energie(System& S,unsigned int N_MC);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int N_spin(3);
	unsigned int N_m(8);
	unsigned int N_MC(1e4);
	P.set("N_spin",N_spin);	
	P.set("N_m",N_m);	
	P.set("N_MC",N_MC);	

	System S(N_spin, N_m, 1);

	Chrono t;
	t.tic();
	std::cout<<energie(S,N_MC)<<std::endl;
	t.tac();
	std::cerr<<t<<" seconde(s)"<<std::endl;
}

double energie(System& S,unsigned int N_MC){
	double ratio(0.0), energie(0.0);
	unsigned int i(0);
	//Save steps("analyse-2d.dat");
	while(i<N_MC){
		S.swap();
		ratio = S.compute_ratio();
		ratio *= ratio;
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			i++;
			S.update_state();
			//steps << S << " " << S.det() <<Save::endl;
			for(unsigned int j(0);j<S.N_site;j++){
				for(unsigned int d(0);d<S.dim;d++){
					S.swap(j,S.nts[S.dim*j+d]);
					energie += S.compute_ratio();
				}
			}
		}
	}
	return -energie/(S.N_site * N_MC); // sign(permutation) => ratio det tjrs - ??? 
}


