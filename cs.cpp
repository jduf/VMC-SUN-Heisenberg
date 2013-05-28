/*!  @file cs.cpp */

#include "CreateSystem.hpp"
#include "Parseur.hpp"

#include <string>
#include <vector>

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	unsigned int N_spin(0), N_m(0), N_n(0);
	bool is_complex(false);
	P.set("N_spin",N_spin);	
	P.set("N_m",N_m);	
	P.set("N_n",N_n);
	P.set("complex",is_complex);
	if(P.n_args()==4){
		if(N_n == 2){
			CreateSystem<double> CS_P(N_m,N_spin,N_n,1);
		} else {
			if(is_complex){
				CreateSystem<std::complex<double> > CS_P(N_m,N_spin,N_n,1);
				CreateSystem<std::complex<double> > CS_AP(N_m,N_spin,N_n,-1);
			} else {
				CreateSystem<double> CS_P(N_m,N_spin,N_n,1);
				CreateSystem<double> CS_AP(N_m,N_spin,N_n,-1);
			}
		}
	} else {
		double td;
		std::string hopfile;
		P.set("td",td);
		P.set("hopfile",hopfile);
		if(is_complex){
			std::cerr<<"main bug"<<std::endl;
		} else {
			CreateSystem<double> CS_AP(N_m,N_spin,N_n,td,hopfile);
		}
	}
}

