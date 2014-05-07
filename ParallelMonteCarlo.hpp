#ifndef DEF_PARALLELMONTECARLO
#define DEF_PARALLELMONTECARLO

#include "MonteCarlo.hpp"

template<typename Type>
class ParallelMonteCarlo {
	public:
		ParallelMonteCarlo(CreateSystem* CS, std::string path, unsigned int nthreads, unsigned int tmax, unsigned int type);
		void run();

		double get_energy() const {return E_;}

	private:
		CreateSystem* CS_;
		unsigned int const tmax_;	
		unsigned int const nruns_;
		unsigned int run_;
		Data<double> E_;
		Vector<double> corr_;
		Vector<double> long_range_corr_;
		Write results_file_;

		unsigned int type_;
};

template<typename Type>
ParallelMonteCarlo<Type>::ParallelMonteCarlo(CreateSystem* CS, std::string path, unsigned int nruns, unsigned int tmax, unsigned int type):
	CS_(CS),
	tmax_(tmax),
	nruns_(nruns),
	//corr_(CS_->get_num_links(),0.0),
	//long_range_corr_(CS->get_n()/3-1,0),
	results_file_(path+CS_->get_filename()+".jdbin"),
	type_(type)
{
	results_file_("type of simulation",type_);
	results_file_("number of simulations runned",nruns_);
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
#pragma omp parallel for 
	for(unsigned int i=0; i<nruns_; i++){
		MonteCarlo<double> sim(CS_,tmax_,type_);
		sim.run();
		E_ += (sim.get_system())->get_energy();
		//corr_ += (sim.get_system())->get_corr();
		//long_range_corr_ += sim.get_long_range_corr();
#pragma omp critical
		{
			(*sim.get_system())>>results_file_;
		}
	}

	RST rst_mean_results;
	rst_mean_results.title("Mean results (status>2)","-");
	results_file_.add_to_header(rst_mean_results.get());
	//results_file_("E (energy per site)",E_);
	//results_file_("DeltaE (absolute error)",DeltaE_);
	//results_file_("corr (correlation on links)",corr_);
	//results_file_("long_range_corr (long range correlation)",long_range_corr_);
}
#endif
