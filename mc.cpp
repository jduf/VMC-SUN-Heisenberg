/*!  @file mc.cpp */

#include "Parseur.hpp"
#include "Read.hpp"
#include "MonteCarlo.hpp"

#include<omp.h>

void run(Parseur& P);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	run(P);
}

void run(Parseur& P){
	std::string filename(P.get<std::string>("sim"));
	if(!P.status()){
		unsigned int nthreads(2);
		P.set("nthreads",nthreads);

		unsigned int N_spin(0), N_m(0);
		bool is_complex(false);

		Read r(filename.c_str());
		r>>is_complex>>N_spin>>N_m;

		Matrix<unsigned int> sts;
		Matrix<double> H;
		r>>sts>>H;
		if(is_complex){
			std::cerr<<"the simulation will be lunched for complex numbers"<<std::endl;
			Matrix<std::complex<double> > T;
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
			Matrix<double> T;
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
		std::cerr<<"main : mc -sim filename.jdbin [-nthreads n]"<<std::endl;
	}
}

