#ifndef DEF_PARALLELMONTECARLO
#define DEF_PARALLELMONTECARLO

#include "MonteCarlo.hpp"

template<typename Type>
class ParallelMonteCarlo {
	public:
		ParallelMonteCarlo(CreateSystem* CS, unsigned int nthreads, unsigned int tmax, unsigned int Nmaxsteps=1e9);
		void run();
		void save(Write& w) const;

		double get_energy() const {return E_;}

	private:
		CreateSystem* CS_;
		unsigned int const tmax_;	
		unsigned int const Nmaxsteps_;	
		unsigned int const nruns_;
		unsigned int run_;
		double E_;
		double DeltaE_;
		Write results_file_;
		Write trash_file_;

		void handle_errors(Vector<double> const& E,Vector<double> const& DeltaE);
};

template<typename Type>
ParallelMonteCarlo<Type>:: ParallelMonteCarlo(CreateSystem* CS, unsigned int nruns, unsigned int tmax, unsigned int Nmaxsteps):
	CS_(CS),
	tmax_(tmax),
	Nmaxsteps_(Nmaxsteps),
	nruns_(nruns),
	E_(0.0),
	DeltaE_(0.0),
	results_file_(CS_->get_filename()+".jdbin"),
	trash_file_(CS_->get_filename()+"-trash.dat")
{
	results_file_("Created by the ParallelMonteCarlo class",0);
	results_file_("nruns (number of simulations runned)",nruns_);
	RST rst_param;
	rst_param.title("Input","-");
	results_file_.add_to_header(rst_param.get());
	CS_->save(results_file_);
	RST rst_results;
	rst_results.title("Results","-");
	results_file_.add_to_header(rst_results.get());
}

template<typename Type>
void ParallelMonteCarlo<Type>::run(){
	Vector<double> E(nruns_);
	Vector<double> DeltaE(nruns_);
#pragma omp parallel for 
	for(unsigned int i=0; i<nruns_; i++){
		MonteCarlo<double> sim(CS_,tmax_,Nmaxsteps_);
		sim.run();
		E(i) = sim.get_energy();
		DeltaE(i) = sim.get_error();
#pragma omp critical
		{
			sim.save(results_file_);
		}
	}
	handle_errors(E,DeltaE);
}

template<typename Type>
void ParallelMonteCarlo<Type>::save(Write& w) const {
	w("E (energy per site)",E_);
	w("DeltaE (absolute error)",DeltaE_);
}

template<typename Type>
void ParallelMonteCarlo<Type>::handle_errors(Vector<double> const& E,Vector<double> const& DeltaE){
	double m(E.mean());
	unsigned int todel(0);
	for(unsigned int i(0);i<E.size();i++){
		if(std::abs(E(i)-m) > 2*DeltaE(i)){ todel++; }
	}
	if(todel>0 && E.size()-todel>2){
		Vector<double> tmpE(E.size()-todel);
		Vector<double> tmpDeltaE(DeltaE.size()-todel);
		Vector<double> trash(todel);
		unsigned int j(0);
		unsigned int k(0);
		for(unsigned int i(0);i<E.size();i++){
			if(std::abs(E(i)-m) <= 2*DeltaE(i)){
				tmpE(j) = E(i);
				tmpDeltaE(j) = DeltaE(i);
				j++;
			} else {
				trash(k) = E(i);
				k++;
			}
		}
		for(unsigned int j(0);j<todel;j++){
			trash_file_<<CS_->get_param()<<" "<<trash(j)<<Write::endl;
		}
		handle_errors(tmpE,tmpDeltaE);
	} else {
		E_ = E.mean(); 
		DeltaE_ = DeltaE.mean()/sqrt(1.0*DeltaE.size());
	}
}
#endif
