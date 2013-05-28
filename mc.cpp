/*!  @file mc.cpp */

#include "Read.hpp"
#include "Matrice.hpp"
#include "Array2D.hpp"
#include "MonteCarlo.hpp"
#include "Parseur.hpp"

#include<string>
#include<iostream>
#include<cmath>
#include<omp.h>

void run(Parseur& P);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	run(P);
}

void run(Parseur& P){
	if(P.n_args()==1 || P.n_args()==2){
		std::string filename("");
		unsigned int nthreads(2);
		P.set("sim",filename);
		P.set("nthreads",nthreads);

		unsigned int N_spin(0), N_m(0), N_n(0);
		bool is_complex(false);

		Read r(filename.c_str());
		r>>is_complex>>N_spin>>N_m>>N_n;

		Array2D<unsigned int> sts(N_spin*N_m*N_n/2,2);
		Matrice<double> H(N_m*N_spin);
		r>>sts>>H;
		if(is_complex){
			std::cerr<<"the simulation will be lunched for complex numbers"<<std::endl;
			Matrice<std::complex<double> > T(N_m*N_spin);
			MonteCarlo<std::complex<double> > sim(filename,nthreads);
			r>>T;
			sim.init(N_spin,N_m,H,sts,T);
			//omp_set_nested(1);
#pragma omp parallel num_threads(nthreads)
			{
				sim.run(omp_get_thread_num());
			}
			sim.save();
		} else {
			std::cerr<<"the simulation will be lunched for real numbers"<<std::endl;
			Matrice<double> T(N_m*N_spin);
			MonteCarlo<double> sim(filename,nthreads);
			r>>T;
			sim.init(N_spin,N_m,H,sts,T);
			//omp_set_nested(1);
#pragma omp parallel num_threads(nthreads)
			{
				sim.run(omp_get_thread_num());
			}
			sim.save();
		}
	} else {
		std::cerr<<"main : mc sim filename.jdbin nthreads n"<<std::endl;
	}
}

