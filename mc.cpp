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

		std::string wf;
		unsigned int n(0), N(0), m(0);
		double bc(0);
		Matrix<unsigned int> sts;

		Read r(filename.c_str());
		r>>wf>>N>>m>>sts;
		if( wf == "chain" ){
			std::cerr<<"1D chain"<<std::endl;
			Matrix<double> EVec;
			MonteCarlo<double> sim(filename,nthreads);
			r>>EVec>>bc;
#pragma omp parallel num_threads(nthreads)
			{
				sim.init(N,m,sts,EVec,omp_get_thread_num());
				sim.run(omp_get_thread_num());
			}
			Write result(filename+".dat");
			result<<"%N n N_samples E_persite Delta_e Status bc"<<Write::endl;
			for(unsigned int thread(0);thread<nthreads;thread++){
				result<<N
					<<" "<<n;
				sim.save_in_file(result,thread);
				result<<" "<<bc
					<<Write::endl;
			}
		}
		if( wf == "fermi" ){
			std::cerr<<"the simulation will be lunched for real numbers"<<std::endl;

			unsigned int N_row(0),N_col(0);

			Matrix<double> EVec;
			MonteCarlo<double> sim(filename,nthreads);
			r>>EVec>>bc>>N_row>>N_col;
#pragma omp parallel num_threads(nthreads)
			{
				sim.init(N,m,sts,EVec,omp_get_thread_num());
				sim.run(omp_get_thread_num());
			}
			Write result(filename+".dat");
			result<<"%N n N_samples E_persite Delta_e Status bc"<<Write::endl;
			for(unsigned int thread(0);thread<nthreads;thread++){
				result<<N
					<<" "<<n;
				sim.save_in_file(result,thread);
				result<<" "<<bc
					<<Write::endl;
			}
		}
		if( wf == "mu" ){
			std::cerr<<"the simulation will be lunched for real numbers"<<std::endl;

			unsigned int N_row(0),N_col(0);
			double mu(0);

			Matrix<double> EVec;
			MonteCarlo<double> sim(filename,nthreads);
			r>>EVec>>bc>>N_row>>N_col>>mu;
#pragma omp parallel num_threads(nthreads)
			{
				sim.init(N,m,sts,EVec,omp_get_thread_num());
				sim.run(omp_get_thread_num());
			}
			Write result(filename+".dat");
			result<<"%N n N_samples E_persite Delta_e Status bc mu"<<Write::endl;
			for(unsigned int thread(0);thread<nthreads;thread++){
				result<<N
					<<" "<<n;
				sim.save_in_file(result,thread);
				result<<" "<<bc
					<<" "<<mu
					<<Write::endl;
			}
		}
		if( wf == "csl" ){
			std::cerr<<"CSL on the square lattice"<<std::endl;

			unsigned int N_row(0),N_col(0);

			Matrix<std::complex<double> > EVec;
			MonteCarlo<std::complex<double> > sim(filename,nthreads);
			r>>EVec>>bc>>N_row>>N_col;
#pragma omp parallel num_threads(nthreads)
			{
				sim.init(N,m,sts,EVec,omp_get_thread_num());
				sim.run(omp_get_thread_num());
			}
			Write result(filename+".dat");
			result<<"%N n N_samples E_persite Delta_e Status bc"<<Write::endl;
			for(unsigned int thread(0);thread<nthreads;thread++){
				result<<N
					<<" "<<n;
				sim.save_in_file(result,thread);
				result<<" "<<bc
					<<Write::endl;
			}
		}
	} else {
		fprintf(stderr,"main : ./mc -sim filename.jdbin -nthreads n\n");
	}
}
