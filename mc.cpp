/*!
@file mc.cpp
*/

#include "Read.hpp"
#include "Matrice.hpp"
#include "Array2D.hpp"
#include "MonteCarlo.hpp"
#include "Parseur.hpp"

#include<string>
#include<iostream>
#include<omp.h>

void run(Parseur& P);
void test(Parseur& P);
void test_real(   std::string filename, double* rnd_sy, unsigned int rnd_sy_size, double* rnd_ra, unsigned int rnd_ra_size, double* rnd_sw, unsigned int rnd_sw_size);
void test_complex(std::string filename, double* rnd_sy, unsigned int rnd_sy_size, double* rnd_ra, unsigned int rnd_ra_size, double* rnd_sw, unsigned int rnd_sw_size);

int main(int argc, char* argv[]){
	Parseur P(argc,argv);
	//std::cerr<<"vérifer abs, fabs, abs(complex)"<<std::endl;
	//std::cerr<<"vérifer det et compute_ration pour les complex"<<std::endl;
	run(P);
	//test(P);
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
			omp_set_nested(1);
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
			omp_set_nested(1);
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

void test(Parseur& P){
	if(P.n_args()==1 || P.n_args()==2){
		std::string filename1("");
		std::string filename2("");
		P.set("simc",filename1);
		P.set("simr",filename2);
		
		unsigned int rnd_sw_size(1e5);
		unsigned int rnd_ra_size(1e2);
		unsigned int rnd_sy_size(100);
		Rand rnd_sw(rnd_sw_size);
		Rand rnd_ra(rnd_ra_size);
		Rand rnd_sy(rnd_sy_size);
		double* static_rnd_sw(new double[rnd_sw_size]);	
		double* static_rnd_ra(new double[rnd_ra_size]);	
		double* static_rnd_sy(new double[rnd_sy_size]);	
		for(unsigned int i(0);i<rnd_sw_size;i++){
			static_rnd_sw[i] = rnd_sw.get();
		}
		for(unsigned int i(0);i<rnd_ra_size;i++){
			static_rnd_ra[i] = rnd_ra.get();
		}
		for(unsigned int i(0);i<rnd_sy_size;i++){
			static_rnd_sy[i] = rnd_sy.get();
		}
		
		test_complex(filename1, static_rnd_sy, rnd_sy_size, static_rnd_ra, rnd_ra_size, static_rnd_sw, rnd_sw_size);
		test_real(   filename2, static_rnd_sy, rnd_sy_size, static_rnd_ra, rnd_ra_size, static_rnd_sw, rnd_sw_size);
	} else {
		std::cerr<<"main : mc sim1 filename-c.jdbin sim2 filename-r.jdbin"<<std::endl;
	}
}

void test_real(std::string filename, double* rnd_sy, unsigned int rnd_sy_size, double* rnd_ra, unsigned int rnd_ra_size, double* rnd_sw, unsigned int rnd_sw_size){
	unsigned int N_spin(0), N_m(0), N_n(0);
	bool is_complex(false);

	Read r(filename.c_str());
	r>>is_complex>>N_spin>>N_m>>N_n;

	Array2D<unsigned int> sts(N_spin*N_m*N_n/2,2);
	Matrice<double> H(N_m*N_spin);
	r>>sts>>H;
	if(is_complex){
		std::cerr<<"should be a real matrix"<<std::endl;
	} else {
		Matrice<double> T(N_m*N_spin);
		MonteCarlo<double> sim(filename,1);
		r>>T;
		sim.init_test(N_spin,N_m,H,sts,T,rnd_sy,rnd_sy_size);
		sim.run_test(rnd_sw,rnd_sw_size,rnd_ra,rnd_ra_size);
		//sim.save_test();
	}
}

void test_complex(std::string filename, double* rnd_sy, unsigned int rnd_sy_size, double* rnd_ra, unsigned int rnd_ra_size, double* rnd_sw, unsigned int rnd_sw_size){
	unsigned int N_spin(0), N_m(0), N_n(0);
	bool is_complex(true);

	Read r(filename.c_str());
	r>>is_complex>>N_spin>>N_m>>N_n;

	Array2D<unsigned int> sts(N_spin*N_m*N_n/2,2);
	Matrice<double> H(N_m*N_spin);
	r>>sts>>H;
	if(is_complex){
		Matrice<std::complex<double> > T(N_m*N_spin);
		MonteCarlo<std::complex<double> > sim(filename,1);
		r>>T;
		sim.init_test(N_spin,N_m,H,sts,T,rnd_sy,rnd_sy_size);
		sim.run_test(rnd_sw,rnd_sw_size,rnd_ra,rnd_ra_size);
		//sim.save_test();
	} else {
		std::cerr<<"should be a complex matrix"<<std::endl;
	}

}
