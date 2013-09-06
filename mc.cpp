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
		unsigned int nthreads(1);
		P.set("nthreads",nthreads);

		unsigned int N_site(0), N_spin(0), N_m(0),N_row(0),N_col(0);
		bool is_complex(false);
		double bc(0);
		double mu(0);
		Matrix<unsigned int> sts;

		Read r(filename.c_str());
		r>>is_complex>>N_spin>>N_m>>sts;
		N_site = N_m*N_spin;

		if(is_complex){
			std::cerr<<"the simulation will be lunched for complex numbers"<<std::endl;
			Matrix<std::complex<double> > EVec;
			MonteCarlo<std::complex<double> > sim(filename,nthreads);
			r>>EVec>>bc>>N_row>>N_col>>mu;
#pragma omp parallel num_threads(nthreads)
			{
				sim.init(N_spin,N_m,sts,EVec,omp_get_thread_num());
				sim.run(omp_get_thread_num());
			}
			Write result(filename+".dat");
			result<<"%N_spin N_site bc N_samples E_persite Delta_e Status"<<Write::endl;
			for(unsigned int thread(0);thread<nthreads;thread++){
				result<<N_spin
					<<" "<<N_site
					<<" "<<bc
					<<" "<<mu;
				sim.save_in_file(result,thread);
				result<<Write::endl;
			}
		} else {
			std::cerr<<"the simulation will be lunched for real numbers"<<std::endl;
			Matrix<double> EVec;
			MonteCarlo<double> sim(filename,nthreads);
			r>>EVec>>bc>>N_row>>N_col>>mu;
#pragma omp parallel num_threads(nthreads)
			{
				sim.init(N_spin,N_m,sts,EVec,omp_get_thread_num());
				sim.run(omp_get_thread_num());
			}
			Write result(filename+".dat");
			result<<"%N_spin N_site bc N_samples E_persite Delta_e Status"<<Write::endl;
			for(unsigned int thread(0);thread<nthreads;thread++){
				result<<N_spin
					<<" "<<N_site
					<<" "<<bc
					<<" "<<mu;
				sim.save_in_file(result,thread);
				result<<Write::endl;
			}
		}
	} else {
		fprintf(stderr,"main : ./mc -sim filename.jdbin -nthreads n\n");
	}
}
