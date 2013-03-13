/*!
@file mc.cpp
*/

#include "Chrono.hpp"
#include "Read.hpp"
#include "System.hpp"
#include "Matrice.hpp"
#include "Array2D.hpp"
#include "MonteCarlo.hpp"

#include<string>
#include<iostream>
#include<cstdlib>

int main(int argc, char* argv[]){
	//std::cerr<<"vérifer abs, fabs, abs(complex)"<<std::endl;
	//std::cerr<<"vérifer det et compute_ration pour les complex"<<std::endl;
	if(argc==2){
		Chrono t;
		t.tic();
		std::string filename(argv[1]);

		unsigned int N_spin(0), N_m(0), N_n(0);
		bool is_complex(false);

		Read r(filename.c_str());
		r>>is_complex>>N_spin>>N_m>>N_n;

		Array2D<unsigned int> sts(N_spin*N_m*N_n/2,2);
		Matrice<double> H(N_m*N_spin);
		r>>sts>>H;
		if(is_complex){
			Matrice<std::complex<double> > T(N_m*N_spin);
			r>>T;
			System<std::complex<double> > S(N_spin,N_m,T);
			MonteCarlo<std::complex<double> > sim(&S,H,sts);
			sim.run();
			t.tac();
			std::cout<<" "<<t<<std::endl;
		} else {
			Matrice<double> T(N_m*N_spin);
			r>>T;
			System<double> S(N_spin,N_m,T);
			MonteCarlo<double> sim(&S,H,sts);
			sim.run();
			t.tac();
			std::cout<<" "<<t<<std::endl;
		}
	} else {
		std::cerr<<"main : need one argument -> filename.jdbin"<<std::endl;
	}
}
