#ifndef DEF_PARALLELMONTECARLO
#define DEF_PARALLELMONTECARLO

#include "MonteCarlo.hpp"

template<typename Type>
class ParallelMonteCarlo {
	public:
		ParallelMonteCarlo(MCSystem<Type>* S, unsigned int nthreads, unsigned int tmax);
		/*!Run all the Monte-Carlo algorithm and save each one in file w*/
		void run(IOFiles& w);
		/*!Save the averaged mesures in w*/
		void save(IOFiles& w);

	private:
		MCSystem<Type>* S_;
		unsigned int const tmax_;	
		unsigned int const nruns_;
		Data<double> E_;
		DataSet<double> corr_;
		DataSet<double> long_range_corr_;
};

template<typename Type>
ParallelMonteCarlo<Type>::ParallelMonteCarlo(MCSystem<Type>* S, unsigned int nruns, unsigned int tmax):
	S_(S),
	tmax_(tmax),
	nruns_(nruns)
{
	E_.set_conv(true);
	//corr_.set(GS->get_n(),true);
	//if(type == 2){
		//long_range_corr_.set(GS->get_n()/3,true);
	//}
}

template<typename Type>
void ParallelMonteCarlo<Type>::run(IOFiles& w){

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
