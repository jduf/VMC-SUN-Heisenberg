#ifndef DEF_PARALLELMONTECARLO
#define DEF_PARALLELMONTECARLO

#include "MonteCarlo.hpp"

template<typename Type>
class ParallelMonteCarlo {
	public:
		ParallelMonteCarlo(CreateSystem* CS, unsigned int nthreads, unsigned int tmax, unsigned int type);
		/*!Run all the Monte-Carlo algorithm and save each one in file w*/
		void run(IOFiles& w);
		/*!Save the averaged mesures in w*/
		void save(IOFiles& w);

	private:
		CreateSystem* CS_;
		unsigned int const tmax_;	
		unsigned int const nruns_;
		unsigned int const type_;
		Data<double> E_;
		DataSet<double> corr_;
		DataSet<double> long_range_corr_;
};

template<typename Type>
ParallelMonteCarlo<Type>::ParallelMonteCarlo(CreateSystem* CS, unsigned int nruns, unsigned int tmax, unsigned int type):
	CS_(CS),
	tmax_(tmax),
	nruns_(nruns),
	type_(type)
{
	E_.set_conv(true);
	corr_.set(CS->get_n(),true);
	if(type == 2){
		long_range_corr_.set(CS->get_n()/3,true);
	}
}

template<typename Type>
void ParallelMonteCarlo<Type>::run(IOFiles& w){
#pragma omp parallel for 
	for(unsigned int i=0; i<nruns_; i++){
		MonteCarlo<double> sim(CS_,tmax_,type_);
		sim.run();
#pragma omp critical
		{
			E_.add_sample((sim.get_system())->get_energy());
			corr_.add_sample((sim.get_system())->get_corr());
			long_range_corr_.add_sample((sim.get_system())->get_long_range_corr());
			(sim.get_system())->save(w);
		}
	}

	E_.complete_analysis();
	corr_.complete_analysis();
	long_range_corr_.complete_analysis();
}

template<typename Type>
void ParallelMonteCarlo<Type>::save(IOFiles& w){
	RST rst_mean_results;
	rst_mean_results.title("Mean results (status>2)","-");
	w.add_to_header(rst_mean_results.get());
	w("E (energy per site)",E_);
	w("corr (correlation on links)",corr_);
	w("lon_range_corr (long range correlation)",long_range_corr_);
}
#endif
