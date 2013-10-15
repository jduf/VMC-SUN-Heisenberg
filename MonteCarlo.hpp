#ifndef DEF_MONTECARLO
#define DEF_MONTECARLO

#include "System.hpp"
#include "Write.hpp"
#include "Chrono.hpp"

#include <vector>

//{Description
/*! Class MonteCarlo
 *
 * Is created simply by giving the name of the output file and the number of
 * independant thread that sould be lunched. An independant System will be
 * created on each thread and therefore, each System should be initialized
 * using the init method.  For each thread, the class let System evolve using
 * its update(), swap() and * ratio() methods.  Thus System.hpp need to have
 * those methods
 */
//}
template <typename Type>
class MonteCarlo{
	public:
		/*!Allocate memory with new[] and starts the chronometer*/
		MonteCarlo(std::string filename, unsigned int const& nthreads); 
		/*!Free allocated memory with delete[] and stops the chronometer*/
		~MonteCarlo();

		//{Public methods, core of the Monte-Carlo simulation
		/*!Lunches an independent simulation on a thread. When one of the
		 * threads is stopped by the test_convergence() method, all thread will
		 * be stopped at the same time and the results for all threads are
		 * computed
		 * */
		//}
		void run(unsigned int const& thread);
		/*!Initializes a different System for each thread*/
		void init(unsigned int const& N_spin, unsigned int const& N_m, Matrix<unsigned int> const& sts, Matrix<Type> const& EVec, unsigned int thread);
		/*!Saves the essential data in the "result" file*/
		void save_in_file(Write& w, unsigned int thread);

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
		void test_convergence(unsigned int const& thread);
		/*!Proceed to a binning analysis and, for each bin, it save the variance */
		void binning(std::vector<double>& d, unsigned int const& thread);
		/*!Computes the mean of a std::vector<double> */
		double mean(std::vector<double> const& v);
		/*!Compute the variance of a std::vector<double> */
		double delta(std::vector<double> const& v, double const& m);
	
		unsigned int const nthreads; //!< Number of independant chains that are lunched
		unsigned int const N_MC; //!< Number of measure to do before doing a binning analysis
		unsigned int const time_limit; //!< Time limit in second
		bool keep_measuring; //!< True if the code runs
		Chrono stop; //!< To stop the simulation after time_limit seconds
		Write output; //!< Text file of name "filename-MC.out" that stores all the output of the MC algorithm
		System<Type>* S; //!< Pointer to a system 
		double* E; //!< Value that the MC algorithm tries to compute
		double* err; //!< Error on E
		std::vector<double>* sampling; //!< Stores all the values that MC considers
		Matrix<unsigned int>* lattice;
		unsigned int* status; //!< Not Lunched:0 Lunched:1 Successful:2 Time elapsed:3
};

/*constructors and destructor*/
/*{*/
template<typename Type>
MonteCarlo<Type>::MonteCarlo(std::string filename, unsigned int const& nthreads):
	nthreads(nthreads),
	N_MC(1e4),
	time_limit(nthreads*5*60),
	keep_measuring(true),
	output(filename+"-MC.out"),
	S(new System<Type>[nthreads]),
	E(new double[nthreads]),
	err(new double[nthreads]),
	sampling(new std::vector<double>[nthreads]),
	lattice(new Matrix<unsigned int>[nthreads]),
	status(new unsigned int[nthreads])
{
	stop.tic();
}

template<typename Type>
MonteCarlo<Type>::~MonteCarlo(){
	delete[] S;
	delete[] lattice;
	delete[] sampling;
	delete[] E;
	delete[] err;
	stop.tac();
}
/*}*/

/*public methods*/
/*{*/
template<typename Type>
void MonteCarlo<Type>::init(unsigned int const& N_spin, unsigned int const& N_m, Matrix<unsigned int> const& sts, Matrix<Type> const& EVec, unsigned int thread) {
	status[thread] = S[thread].init(N_spin,N_m,sts,EVec,thread);
	if(status[thread]){
		std::cerr<<"thermalization "<<std::flush;
		unsigned int i(0);
		while( i<1e5 ){
			S[thread].swap();
			if( std::abs(S[thread].ratio())>1e-10 ){
				i++;
				S[thread].update();
			}
		}
		std::cerr<<"completed"<<std::endl;
		lattice[thread].set(N_spin*N_m,N_spin,0);
	}
}

template<typename Type>
void MonteCarlo<Type>::run(unsigned int const& thread){
	if(status[thread]){
		unsigned int i(0);
		double E_config(0);
		double ratio(0.0);
		Rand rnd(1e4,thread);
		bool bin(true);
		if(bin){
			//do{
				//S[thread].swap();
				//ratio = norm_squared(S[thread].ratio());
				//if( ratio > 1.0 || rnd.get() < ratio ){
					//S[thread].update();
					//S[thread].measure(E_config,lattice[thread]);
					//E[thread] += E_config;
					//sampling[thread].push_back(E_config);
					//i++;
					//if(i >= N_MC) {
						//i=0;
						//test_convergence(thread);
						//output<<sampling[thread].size()
							//<<" "<<E[thread] / (sampling[thread].size() * S[thread].N_site)
							//<<" "<<err[thread]
							//<<Write::endl;
					//}
				//}
			//} while(keep_measuring);
			do{
				S[thread].swap();
				ratio = norm_squared(S[thread].ratio());
				if( ratio > 1.0 || rnd.get() < ratio ){
					S[thread].update();
					S[thread].measure(E_config,lattice[thread]);
					E[thread] += E_config;
					sampling[thread].push_back(E_config);
				} else {
					E[thread] += E_config;
					sampling[thread].push_back(E_config);
				}
				i++;
				if(i >= N_MC) {
					i=0;
					test_convergence(thread);
					output<<sampling[thread].size()
						<<" "<<E[thread] / (sampling[thread].size() * S[thread].N_site)
						<<" "<<err[thread]
						<<Write::endl;
				}
			} while(keep_measuring);
		} else {
			double a(0.0);
			do{
				S[thread].swap();
				ratio = norm_squared(S[thread].ratio());
				a = rnd.get();
				i++;
				if(ratio>1 || a <ratio){
					S[thread].update();
					S[thread].measure(E_config,lattice[thread]);
					E[thread] += E_config;
					sampling[thread].push_back(E_config);
					std::cout<<i<<" "<<a<<" "<<ratio<<" "<<E_config<<std::endl;
				} else {
					E[thread] += E_config;
				}
			} while(i<1e4); 
			output<<lattice[thread]<<Write::endl;
		}
		test_convergence(thread);
	}
}

template<typename Type>
void MonteCarlo<Type>::save_in_file(Write& w, unsigned int thread){
	w<<" "<<sampling[thread].size()
		<<" "<<E[thread] / (sampling[thread].size() * S[thread].N_site)
		<<" "<<err[thread]
		<<" "<<status[thread];
}
/*}*/

/*private methods*/
/*{*/
template<typename Type>
void MonteCarlo<Type>::test_convergence(unsigned int const& thread){
	if( std::abs( E[thread] / sampling[thread].size() )  > 1e3 ){ 
		E[thread] = 0.0;
		sampling[thread].clear();
		std::cerr<<"the simulation is restarted, bad initial condition"<<std::endl;
	}  else {
		if(keep_measuring && stop.time_limit_reached(time_limit)){
			keep_measuring = false;
			status[thread] = 3;
			std::cerr<<"the simulation was stopped because it reached the time limit"<<std::endl;
		}

		std::vector<double> d(0);
		binning(d,thread);

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
		err[thread] = d[d.size()-1] / S[thread].N_site; 
		if(cond/m<1e-3){ 
			keep_measuring = false; 
			status[thread] = 2;
		}
	}
}

template<typename Type>
void MonteCarlo<Type>::binning(std::vector<double>& d, unsigned int const& thread){
	std::vector<double> bin(sampling[thread]);
	std::vector<double> bin2(bin.size()/2);
	unsigned int l(0);
	while(pow(2,l+1)*100<sampling[thread].size()){
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
