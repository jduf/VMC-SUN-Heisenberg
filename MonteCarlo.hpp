#ifndef DEF_MONTECARLO
#define DEF_MONTECARLO

#include "System.hpp"
#include "Array2D.hpp"
#include "Write.hpp"
#include "Rand.hpp"
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
template <typename T>
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
		void init(unsigned int const& N_spin, unsigned int const& N_m, Matrice<double> const& H, Array2D<unsigned int> const& sts, Matrice<T> EVec);
		/*!Saves the essential data in the "result" file*/
		void save();

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
		double delta(std::vector<double> const& v, double m);

	
		unsigned int const nthreads; //!< Number of independant chains that are lunched
		unsigned int const time_limit; //!< Time limit in second
		unsigned int const N_MC; //!< Number of measure to do before doing a binning analysis
		System<T>* S; //!< Pointer to a system 
		double* E; //!< Value that the MC algorithm tries to compute
		double* err; //!< Error on E
		unsigned int status; //!< Successful:0 Not lunched:1  Time elapsed:2
		std::vector<double>* sampling; //!< Stores all the values that MC considers
		std::string filename; //!< Name of the output file
		Write output; //!< Text file of name "filename.out" that stores all the output of the MC algorithm
		Write result; //!< Text file of name "filename.dat" that stores the final result
		bool keep_measuring; //!< True if the code runs
		Chrono stop; //!< To stop the simulation after time_limit seconds
};

/*constructors and destructor*/
/*{*/
template<typename T>
MonteCarlo<T>::MonteCarlo(std::string filename, unsigned int const& nthreads):
	nthreads(nthreads),
	time_limit(nthreads*3600*2*24),
	N_MC(1e4),
	S(new System<T>[nthreads]),
	E(new double[nthreads]),
	err(new double[nthreads]),
	status(1),
	sampling(new std::vector<double>[nthreads]),
	filename(filename),
	output(filename+".out"),
	result(filename+".dat"),
	keep_measuring(true)
{
	result<<"%N_spin "<<"N_m "<<"N_samples "<<"E_persite "<<"Delta_E "<<"Boundary_condition "<<"Status"<<Write::endl;
	stop.tic();
}

template<typename T>
MonteCarlo<T>::~MonteCarlo(){
	delete[] S;
	delete[] sampling;
	delete[] E;
	delete[] err;
	stop.tac();
}
/*}*/

/*public void methods*/
/*{*/
template<typename T>
void MonteCarlo<T>::run(unsigned int const& thread){
	unsigned int i(1);
	double E_config(0);
	double ratio(0.0);
	Rand rnd(1e4,thread);
	do{
		S[thread].swap();
		ratio = norm_squared(S[thread].ratio());
		if(ratio>1 || rnd.get() <ratio){
			S[thread].update();
			E_config = S[thread].compute_energy();
			E[thread] += E_config;
			sampling[thread].push_back(E_config);
			if(i < N_MC){ i++;}
			else {
				i = 1;
				test_convergence(thread);
			}
		}
	} while(keep_measuring);
#pragma omp critical
	{
		test_convergence(thread);
	}
}

template<typename T>
void MonteCarlo<T>::init(unsigned int const& N_spin, unsigned int const& N_m, Matrice<double> const& H, Array2D<unsigned int> const& sts, Matrice<T> EVec){
	for(unsigned int i(0);i<nthreads;i++){
		S[i].init(N_spin,N_m,H,sts,EVec,i);
	}
}

template<typename T>
void MonteCarlo<T>::save(){
	int parity(0);
	if(filename.find("P") != std::string::npos ){ parity = -1;}
	if(filename.find("A") != std::string::npos ){ parity = 1;}
	for(unsigned int thread(0); thread<nthreads; thread++){
		result<<S[thread].N_spin
			<<" "<<S[thread].N_m
			<<" "<<sampling[thread].size()
			<<" "<<E[thread] / (sampling[thread].size() * S[thread].N_site)
			<<" "<<err[thread]
			<<" "<<parity
			<<" "<<status
			<<Write::endl;
	}
}
/*}*/

/*private void methods*/
/*{*/
template<typename T>
void MonteCarlo<T>::test_convergence(unsigned int const& thread){
	if( fabs( E[thread] / sampling[thread].size() )  > 1e3 ){ 
		std::cerr<<filename<< " : initial condition lead to a wrong value, restarting the simulation (E="<<E[thread]<<")"<<std::endl;
		E[thread] = 0;
		sampling[thread].clear();
	}  else {
		if(keep_measuring && stop.time_limit_reached(time_limit)){
			keep_measuring = false;
			status = 2;
			result<<"% the simulation was stopped because it reached the time limit"<<Write::endl;
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
		if(cond/m<1e-2){ 
			keep_measuring = false; 
			status = 0;
		}
	
		output<<S[thread].N_spin
			<<" "<<S[thread].N_m
			<<" "<<sampling[thread].size()
			<<" "<<E[thread] / (sampling[thread].size() * S[thread].N_site)
			<<" "<<err[thread];
		for(unsigned int i(0);i<d.size();i++){
			output<<" "<<d[i];
		}
		output<<Write::endl;
	}
}

template<typename T>
void MonteCarlo<T>::binning(std::vector<double>& d, unsigned int const& thread){
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
/*}*/

/*private methods with return*/
/*{*/
template<typename T>
double MonteCarlo<T>::mean(std::vector<double> const& v){
	double m(0.0);
	unsigned int N(v.size());
	for(unsigned int i(0);i<N;i++){
		m += v[i];
	}
	return m/N;
}

template<typename T>
double MonteCarlo<T>::delta(std::vector<double> const& v, double m){
	double d(0.0);
	unsigned int N(v.size());
	for(unsigned int i(0);i<N;i++){
		d += (v[i]-m)*(v[i]-m);
	}
	return sqrt(d / (N*(N-1))); 
}
/*}*/

/*double norm_squared(T)*/
/*{*/
template<typename T>
double norm_squared(T x);

template<>
inline double norm_squared(double x){
	return x*x;
}

template<>
inline double norm_squared(std::complex<double> x){
	return std::norm(x);
}
/*}*/
#endif
