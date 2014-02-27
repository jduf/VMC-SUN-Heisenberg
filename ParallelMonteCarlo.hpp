#ifndef DEF_PARALLELMONTECARLO
#define DEF_PARALLELMONTECARLO

#include "MonteCarlo.hpp"

template<typename Type>
class ParallelMonteCarlo {
	public:
		ParallelMonteCarlo(CreateSystem const& CS, double param, unsigned int nthreads, unsigned int tmax);
		void run(unsigned int N_MC=0);
		void save(Write& w);

		double get_energy() const {return E_;}

	private:
		CreateSystem CS_;
		unsigned int tmax_;	
		unsigned int nruns_;
		unsigned int run_;
		double E_;
		double DeltaE_;
		Write results_file_;
};

template<typename Type>
ParallelMonteCarlo<Type>:: ParallelMonteCarlo(CreateSystem const& CS, double param, unsigned int nruns, unsigned int tmax):
	CS_(CS,param),
	tmax_(tmax),
	nruns_(nruns),
	E_(0.0),
	DeltaE_(0.0),
	results_file_(CS_.get_filename()+".jdbin")
{
	results_file_("Created by the ParallelMonteCarlo class",0);
	results_file_("nruns (number of simulations runned)",nruns_);
	RST rst_param;
	rst_param.title("Input","-");
	results_file_.add_to_header(rst_param.get());
	CS_.save(results_file_);
	RST rst_results;
	rst_results.title("Results","-");
	results_file_.add_to_header(rst_results.get());
}

template<typename Type>
void ParallelMonteCarlo<Type>::run(unsigned int N_MC){
	#pragma omp parallel for 
	for(unsigned int i=0; i<nruns_; i++){
		MonteCarlo<double> sim(CS_,tmax_);
		sim.run(N_MC);
		E_ += sim.get_energy();
		DeltaE_ += sim.get_error();
#pragma omp critical
		{
			sim.save(results_file_);
		}
	}
	E_/=nruns_;
	DeltaE_/=nruns_*sqrt(1.0*nruns_);
}

template<typename Type>
void ParallelMonteCarlo<Type>::save(Write& w){
	w("E (energy per site)",E_);
	w("DeltaE (absolute error)",DeltaE_);
}
#endif
