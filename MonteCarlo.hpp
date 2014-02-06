#ifndef DEF_MONTECARLO
#define DEF_MONTECARLO

#include "SystemBosonic.hpp"
#include "SystemFermionic.hpp"
#include "Write.hpp"

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
 *
 * Each time that a class is instanciated, a random number generator is created
 * according to the thread on which the code is running. The same thread number
 * is then transmitted to create the (System*). Once the (*System) is created,
 * it is thermilized.
 *
 * The run(what,N_MC) method lunches the Monte-Carlo simulation.
 */
//}

template <typename Type>
class MonteCarlo{
	public:
		/*!The container must contain a filename to store the executing data,
		 * the parameter required for (*System) and an optional t_max option*/
		MonteCarlo(Container const& input); 
		~MonteCarlo();

		/*!Run the Monte-Carlo algorithme */
		void run(unsigned int what, unsigned int N_MC, Matrix<unsigned int>* lattice = NULL);
		/*!Saves the essential data in the "result" file*/
		Container save();

		/*!Get the energy per site*/
		double get_energy() const { return this->E_/(sampling.size() * S->n_); };
		unsigned int get_status() const { return status; };

	private:
		/*!Forbids the copy constructor*/
		MonteCarlo(MonteCarlo const& mc); 
		/*!Forbids the assignment operator*/
		MonteCarlo const& operator=(MonteCarlo const& mc);

		/*!Find the next configuration and measure it*/
		void next_step(double& E_step);
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
		/*!Proceed to a binning analysis and, for each bin, it saves the variance */
		void binning(std::vector<double>& d);
		/*!Computes the mean of a std::vector<double> */
		double mean(std::vector<double> const& v);
		/*!Compute the variance of a std::vector<double> */
		double delta(std::vector<double> const& v, double const& m);

		System<Type>* S; 		//!< Pointer to a Fermionic or Bosonic System 
		Rand* rnd; 				//!< Pointer to a random number generator
		double E_; 				//!< Value that the MC algorithm tries to compute
		double err; 			//!< Error on E_
		std::vector<double> sampling; //!< Stores all the values that MC considers
		Time stop; 				//!< To stop the simulation after time_limit seconds
		unsigned int t_max;		//!< Time limit in second, by default 5min
		unsigned int status; 	//!< Not Lunched:0 Lunched:1 Successful:2 Time elapsed:3
		bool keep_measuring; 	//!< True if the code runs
		std::string filename_;	//!< Filename for the output
};

/*constructors and destructor*/
/*{*/
template<typename Type>
MonteCarlo<Type>::MonteCarlo(Container const& input):
	S(NULL),
	rnd(NULL),
	E_(0),
	err(0),
	t_max(5*60),
	status(0),
	keep_measuring(true),
	filename_("")
{
	if(input.check("filename")){
		input.get("filename",filename_);
	}
	unsigned int thread(omp_get_thread_num());
	rnd=new Rand(1e4,thread);
	if(input.get<Vector<unsigned int> >("ref")(0)){ S=new SystemFermionic<Type>; }
	else { S=new SystemBosonic<Type>; }
	input.get("t_max",t_max);
	status = S->init(input,thread);
	//std::cout<<"thermalization"<<std::endl;
	if(status){
		double ratio(0.0);
		for(unsigned int i(0);i<1e5;i++){
			S->swap();
			ratio = norm_squared(S->ratio());
			if( ratio > rnd->get() ){
				S->update();
			}
		}
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
void MonteCarlo<Type>::run(unsigned int what, unsigned int N_MC, Matrix<unsigned int>* lattice){
	if(status){
		double E_step(0.0);
		S->measure(E_step);
		switch(what){
			case 1: 
				{
					Write output(filename_+".out");
					do{
						for(unsigned int i(0);i<N_MC;i++){ next_step(E_step); }
						test_convergence();
						output<<sampling.size()
							<<" "<<E_ / (sampling.size() * S->n_)
							<<" "<<err
							<<Write::endl;
					} while(keep_measuring);
				} break;
			case 2:
				{
					for(unsigned int i(0);i<N_MC;i++){ next_step(E_step); }
					test_convergence();
				} break;
			case 3:
				{
					lattice->set(S->n_,S->N_,0); 
					Matrix<unsigned int> corr(S->n_,6,0); 
					for(unsigned int i(0);i<N_MC;i++){ 
						next_step(E_step); 
						/*to check the color organization*/
						for(unsigned int i(0);i<S->n_;i++){ 
							(*lattice)(i,S->s_(i,0))++;
						}
						S->correlation(&corr);
					}
					test_convergence();
					std::cout<<corr<<std::endl;
				} break;
		}
	}
}

template<typename Type>
Container MonteCarlo<Type>::save(){
	Container out(true);
	out.set("N_step",sampling.size());
	out.set("E",E_ / (sampling.size() * S->n_));
	out.set("dE",err);
	out.set("status",status);
	return out;
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
void MonteCarlo<Type>::next_step(double& E_step){
	S->swap();
	if( norm_squared(S->ratio()) > rnd->get() ){
		S->update();
		S->measure(E_step);
	}
	E_ += E_step;
	sampling.push_back(E_step);
}

template<typename Type>
void MonteCarlo<Type>::test_convergence(){
	if( std::abs( E_ / sampling.size() )  > 1e3 ){ 
		E_ = 0.0;
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
