#ifndef DEF_MONTECARLO
#define DEF_MONTECARLO

#include "SystemBosonic.hpp"
#include "SystemFermionic.hpp"
#include "Write.hpp"
#include "Time.hpp"

#include <vector>

//{Description
/*! Class MonteCarlo
 *
 * Implement the Monte-Carlo algorithm. This alorithm let a system (System*)
 * evolves according to a Markov process. To work properly, System* has to be
 * written with :
 * 
 * + System->ratio() : compute the probability to accept the next configuration
 * + System->update() : update the old cufiguration to the new one
 * + System->measure() : measure an observable according to the current
 * configuration
 */
//}
template <typename Type>
class MonteCarlo{
	public:
		/*!*/
		MonteCarlo(Container const& input); 
		~MonteCarlo();

		/*!Run the Monte-Carlo algorithme */
		void run();
		/*!Saves the essential data in the "result" file*/
		Container save();

	private:
		/*!Forbids the copy constructor*/
		MonteCarlo(MonteCarlo const& mc);
		/*!Forbids the assignment operator*/
		MonteCarlo const& operator=(MonteCarlo const& mc);


		//{Private methods that give a shutoff condition
		/*!Stops the simulation when
		 *
		 * - convergence is reached
		 * - time limit is up
		 * - if the energy diverges
		 *
		 * ands save each step in the "output" file
		 */
		//}
		void test_convergence();
		/*!Proceed to a binning analysis and, for each bin, it save the variance */
		void binning(std::vector<double>& d);
		/*!Computes the mean of a std::vector<double> */
		double mean(std::vector<double> const& v);
		/*!Compute the variance of a std::vector<double> */
		double delta(std::vector<double> const& v, double const& m);
	
		System<Type> *S; 		//!< Pointer to a Fermionic or Bosonic System 
		Rand* rnd; 				//!< Pointer to a random number generator
		double E; 				//!< Value that the MC algorithm tries to compute
		double err; 			//!< Error on E
		std::vector<double> sampling; //!< Stores all the values that MC considers
		Time stop; 				//!< To stop the simulation after time_limit seconds
		unsigned int t_max;		//!< Time limit in second, by default 5min
		unsigned int status; 	//!< Not Lunched:0 Lunched:1 Successful:2 Time elapsed:3
		unsigned int const N_MC;//!< Number of measure to do before doing a binning analysis
		bool keep_measuring; 	//!< True if the code runs
		Write output; 			//!< Text file of name "filename-MC.out" that stores all the output of the MC algorithm
};

/*constructors and destructor*/
/*{*/
template<typename Type>
MonteCarlo<Type>::MonteCarlo(Container const& input):
	S(NULL),
	rnd(NULL),
	E(0),
	err(0),
	t_max(5*60),
	status(0),
	N_MC(1e4),
	keep_measuring(true),
	output(input.get<std::string>("filename")+".out")
{
	//unsigned int thread(omp_get_thread_num());
	unsigned int thread(0);
	std::cerr<<"Be sure that each Monte Carlo has its own random number generator, and that this generator is different than the System one"<<std::endl;
	rnd=new Rand(1e4,thread);
	if(input.get<bool>("fermionic")){ S=new SystemFermionic<Type>; }
	else { S=new SystemBosonic<Type>;}
	input.get("t_max",t_max);
	status = S->init(input,thread);
	//if(status){
		//S->print();
		//for(unsigned int i(0);i<2;i++){
			//S->swap();
			//S->ratio();
			//S->update();
			//S->print();
		//}
	//}
	
	if(status){
		std::cerr<<"thermalization (has been changed)"<<std::flush;
		double ratio(0.0);
		for(unsigned int i(0);i<1e5;i++){
			S->swap();
			ratio = norm_squared(S->ratio());
			if( ratio > rnd->get() ){
				S->update();
			}
		}
		std::cerr<<"completed"<<std::endl;
	}
}

template<typename Type>
MonteCarlo<Type>::~MonteCarlo(){
	delete S;
	delete rnd;
}
/*}*/

/*public methods*/
/*{*/
template<typename Type>
void MonteCarlo<Type>::run(){
	if(status){
		unsigned int i(0);
		double E_config(0);
		double ratio(0.0);
		bool bin(false);
		std::cerr<<"there is maybe a problem, for the first iterations, if the new state is rejected, I add a  0 contribution to the energy..."<<std::endl;
		if(bin){
			do{
				S->swap();
				ratio = norm_squared(S->ratio());
				if( ratio > rnd->get() ){
					S->update();
					S->measure(E_config);
				}
				E += E_config;
				sampling.push_back(E_config);
				i++;
				if(i >= N_MC){
					i=0;
					test_convergence();
					output<<sampling.size()
						<<" "<<E / (sampling.size() * S->n_)
						<<" "<<err
						<<Write::endl;
				}
			} while(keep_measuring);
		} else {
			Matrix<unsigned int> lattice(S->n_,S->N_,0); 
			do{
				S->swap();
				ratio = norm_squared(S->ratio());
				if( ratio > rnd->get() ){
					S->update();
					S->measure(E_config);
				}
				E += E_config;
				sampling.push_back(E_config);
				i++;
				/*to check the color organization*/
				for(unsigned int i(0);i<S->n_;i++){ 
					lattice(i,S->s_(i,0))++;
				}
			} while(i<1e5); 
			std::cout<<lattice<<std::endl;
			S->print();
			//for(unsigned int i(0);i<sampling.size();i++){
			//std::cout<<i<<" "<<sampling[i]<<std::endl;
			//}
			test_convergence();
		}
	}
}

template<typename Type>
Container MonteCarlo<Type>::save(){
	Container out(true);
	out.set("N_step",sampling.size());
	out.set("E",E / (sampling.size() * S->n_));
	out.set("dE",err);
	out.set("status",status);
	return out;
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
void MonteCarlo<Type>::test_convergence(){
	if( std::abs( E / sampling.size() )  > 1e3 ){ 
		E = 0.0;
		sampling.clear();
		std::cerr<<"the simulation is restarted, bad initial condition"<<std::endl;
	}  else {
		if(keep_measuring && stop.time_limit_reached(t_max)){
			keep_measuring = false;
			status = 3;
			std::cerr<<"the simulation was stopped because it reached the time limit"<<std::endl;
		}

		std::vector<double> d(0);
		binning(d);

		double v[3];
		double m(0.0);
		double cond(0.0);
		for(unsigned int i(0); i<3; i++){
			v[i] = d[d.size()-1-i];
			m += v[i];
		}
		m /= 3;
		for(unsigned int i(0);i<2;i++){
			cond += (v[i] - m)*(v[i] - m);
		}
		cond += sqrt((cond)/2);
		err = d[d.size()-1] / S->n_; 
		if(cond/m<1e-3){ 
			keep_measuring = false; 
			status = 2;
		}
	}
}

template<typename Type>
void MonteCarlo<Type>::binning(std::vector<double>& d){
	std::vector<double> bin(sampling);
	std::vector<double> bin2(bin.size()/2);
	unsigned int l(0);
	while(pow(2,l+1)*100<sampling.size()){
		d.push_back(delta(bin,mean(bin)));
		l++;
		for(unsigned int i(0);i<bin2.size();i++){
			bin2[i] = (bin[2*i] + bin[2*i+1])/2.0;
		}
		bin.resize(bin2.size());
		bin=bin2;
		bin2.resize(bin2.size()/2);
	}
}

template<typename Type>
double MonteCarlo<Type>::mean(std::vector<double> const& v){
	double m(0.0);
	unsigned int N(v.size());
	for(unsigned int i(0);i<N;i++){
		m += v[i];
	}
	return m/N;
}

template<typename Type>
double MonteCarlo<Type>::delta(std::vector<double> const& v, double const& m){
	double d(0.0);
	unsigned int N(v.size());
	for(unsigned int i(0);i<N;i++){
		d += (v[i]-m)*(v[i]-m);
	}
	return sqrt(d / (N*(N-1))); 
}
/*}*/
#endif
