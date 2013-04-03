/*!
@file cs.cpp
*/

#include "CreateSystem.hpp"
#include "Parseur.hpp"

#include <string>
#include <vector>

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int N_spin(0), N_m(0), N_n(0);
	bool is_complex(false);
	std::string filename;
	P.set("N_spin",N_spin);	
	P.set("N_m",N_m);	
	P.set("N_n",N_n);
	P.set("complex",is_complex);

	if(is_complex){
		CreateSystem<std::complex<double> > CS_P(N_m,N_spin,N_n,-1);
		CreateSystem<std::complex<double> > CS_AP(N_m,N_spin,N_n,1);
	} else {
		CreateSystem<double> CS_P(N_m,N_spin,N_n,-1);
		CreateSystem<double> CS_AP(N_m,N_spin,N_n,1);
	}
}

