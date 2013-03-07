/*!
@file MC.cpp
*/

#include "Chrono.hpp"
#include "Read.hpp"
#include "System.hpp"
#include "Matrice.hpp"
#include "Array2D.hpp"

#include<string>
#include<iostream>
#include<cstdlib>

double energie(System<std::complex<double> >& S,unsigned int N_MC);
double energie(System<double>& S,unsigned int N_MC);

int main(int argc, char* argv[]){
	//std::cerr<<"vérifer abs, fabs, abs(complex)"<<std::endl;
	//std::cerr<<"vérifer det et compute_ration pour les complex"<<std::endl;
	if(argc==2){
		Chrono t;
		t.tic();
		std::string filename(argv[1]);

		unsigned int N_spin(0), N_m(0), N_n(0), N_MC(1e6);
		bool is_complex(false);

		Read r(filename.c_str());
		r>>is_complex>>N_spin>>N_m>>N_n;

		Array2D<unsigned int> sts(N_spin*N_m*N_n/2,2);
		r>>sts;
		if(is_complex){
			Matrice<std::complex<double> > EVec(N_m*N_spin);
			Matrice<std::complex<double> > H(N_m*N_spin);
			r>>H>>EVec;
			System<std::complex<double> > S(N_spin,N_m,sts,H,EVec);
			std::cerr<<is_complex<<" "<<N_n<<" "<<N_spin<<" "<<N_spin*N_m<<" "<<N_MC<<" "<<energie(S,N_MC);
			t.tac();
			std::cerr<<" "<<t<<std::endl;
		} else {
			Matrice<double> EVec(N_m*N_spin);
			Matrice<double> H(N_m*N_spin);
			r>>H>>EVec;
			System<double> S(N_spin,N_m,sts,H,EVec);
			std::cerr<<is_complex<<" "<<N_n<<" "<<N_spin<<" "<<N_spin*N_m<<" "<<N_MC<<" "<<energie(S,N_MC);
			t.tac();
			std::cerr<<" "<<t<<std::endl;
		}
	} else {
		std::cerr<<"main : need one argument -> filename.jdbin"<<std::endl;
	}
}

double energie(System<double>& S,unsigned int N_MC){
	srand(time(NULL)^(getpid()<<16));
	double ratio(0.0);
	unsigned int i(0);
	unsigned int DCT(0);
	double E(0.0);
	double E_state(0.0);
	while(i<N_MC){
		S.swap();
		ratio = S.ratio(false); 
		ratio *= ratio;
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			S.update();
			if(DCT == 0){
				for(unsigned int j(0);j<S.sts.row();j++){
					S.swap(S.sts(j,0),S.sts(j,1));
					E_state += S.ratio(true);
				}
				i++;
				DCT = S.N_site;
				E += E_state;
				std::cout<<E_state<<" "<<E/(i * S.N_site)<<std::endl;
				E_state = 0.0;
			}
			DCT--;
		}
	}
	return E/(S.N_site * N_MC); 
}

double energie(System<std::complex<double> >& S,unsigned int N_MC){
	srand(time(NULL)^(getpid()<<16));
	double ratio(0.0);
	unsigned int i(0);
	unsigned int DCT(0);
	std::complex<double> E(0.0);
	std::complex<double> E_state(0.0);
	while(i<N_MC){
		S.swap();
		ratio = std::norm(S.ratio(false));
		if(ratio>1 || (double)rand()/RAND_MAX <ratio){
			S.update();
			if(DCT == 0){
				for(unsigned int j(0);j<S.sts.row();j++){
					S.swap(S.sts(j,0),S.sts(j,1));
					E_state += S.ratio(true);
				}
				i++;
				DCT = S.N_site;
				E += E_state; 
				std::cout<<std::real(E_state)<<" "<<std::real(E)/(i * S.N_site)<<std::endl;
				E_state = 0.0;
			}
			DCT--;
		}
	}
	return std::real(E)/(S.N_site * N_MC); 
}
