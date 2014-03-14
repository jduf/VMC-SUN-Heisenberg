#ifndef DEF_PARALLELMONTECARLO
#define DEF_PARALLELMONTECARLO

#include "MonteCarlo.hpp"

template<typename Type>
class ParallelMonteCarlo {
	public:
		ParallelMonteCarlo(CreateSystem* CS, std::string path, unsigned int nthreads, unsigned int tmax, unsigned int Nmaxsteps=1e9);
		void run();

		double get_energy() const {return E_;}

	private:
		CreateSystem* CS_;
		unsigned int const tmax_;	
		unsigned int const Nmaxsteps_;	
		unsigned int const nruns_;
		unsigned int run_;
		double E_;
		double DeltaE_;
		Vector<double> corr_;
		Vector<double> long_range_corr_;
		Write results_file_;
};

template<typename Type>
ParallelMonteCarlo<Type>:: ParallelMonteCarlo(CreateSystem* CS, std::string path, unsigned int nruns, unsigned int tmax, unsigned int Nmaxsteps):
	CS_(CS),
	tmax_(tmax),
	Nmaxsteps_(Nmaxsteps),
	nruns_(nruns),
	E_(0.0),
	DeltaE_(0.0),
	corr_(CS_->get_num_links(),0.0),
	long_range_corr_(CS->get_n()/3-1,0),
	results_file_(path+CS_->get_filename()+".jdbin")
{
	results_file_("nruns (number of simulations runned)",nruns_);
	RST rst;
	rst.title("Input","-");
	results_file_.add_to_header(rst.get());
	rst.set();
	CS_->save(results_file_);
	rst.title("Results","-");
	results_file_.add_to_header(rst.get());
}

template<typename Type>
void ParallelMonteCarlo<Type>::run(){
	unsigned int n_ok(0);
#pragma omp parallel for 
	for(unsigned int i=0; i<nruns_; i++){
		MonteCarlo<double> sim(CS_,tmax_,Nmaxsteps_);
		sim.run();
		if(sim.get_status()){
			E_ += sim.get_energy();
			DeltaE_ += sim.get_error();
			corr_ += sim.get_corr();
			long_range_corr_ += sim.get_long_range_corr();
			n_ok++;
		}
#pragma omp critical
		{
			sim.save(results_file_);
		}
	}
	E_ /= n_ok;
	DeltaE_ /= (n_ok*sqrt(n_ok));
	corr_  /= n_ok;
	long_range_corr_  /= n_ok;

	RST rst_mean_results;
	rst_mean_results.title("Mean results (status>2)","-");
	results_file_.add_to_header(rst_mean_results.get());
	results_file_("E (energy per site)",E_);
	results_file_("DeltaE (absolute error)",DeltaE_);
	results_file_("corr (correlation on links)",corr_);
	results_file_("long_range_corr (long range correlation)",long_range_corr_);
}
#endif
