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
 * Takes a pointer to System and let the system evolve using its update(), swap() and
 * ratio() methods. Thus System.hpp need to have those methods
 */
//}
template <typename T>
class MonteCarlo{
	public:
		/*!Constructors */
		MonteCarlo(std::string filename); 
		~MonteCarlo();

		void run(unsigned int const& thread);
		void save(unsigned int const& nthreads);
		void init(unsigned int const& N_spin, unsigned int const& N_m, Matrice<double> const& H, Array2D<unsigned int> const& sts, Matrice<T> EVec, unsigned int const& nthreads);

	private:
		/*!Forbids the copy constructor*/
		MonteCarlo(MonteCarlo const& mc);
		/*!Forbids the assignment operator*/
		MonteCarlo const& operator=(MonteCarlo const& mc);

		void test_convergence(unsigned int const& thread);
		void binning(std::vector<double>& d, unsigned int const& thread);

		/*!Computes the mean of a std::vector<double> */
		double mean(std::vector<double> const& v);

		/*!Compute the variance of a std::vector<double> */
		double delta(std::vector<double> const& v, double m);
	
		System<T>* S; //!< Pointer to a system 
		std::vector<double>* sampling; //<! Stores all the values that MC considers
		double* E; //!< Value that the MC algorithm tries to compute
		double* err; //!< Error on E
		unsigned int const N_MC; //!< Number of measure to do before doing a binning analysis
		Write output; //!< Text file of name "filename.out" that stores all the output of the MC algorithm
		Write result; //!< Text file of name "filename.dat" that stores the final result
		std::string filename; //!< Name of the output file
		bool keep_measuring; //!< True if the code runs
		Chrono stop; //!< To stop the simulation after time_limit seconds
		unsigned int time_limit; //!< Time limit in second
};

template<typename T>
MonteCarlo<T>::MonteCarlo(std::string filename):
	N_MC(1e4),
	output(filename+".out"),
	result(filename+".dat"),
	filename(filename),
	keep_measuring(true),
	time_limit(10)
{
	result<<"%N_spin "<<"N_m "<<"N_samples "<<"E_persite "<<"Delta_E "<<Write::endl;
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

template<typename T>
void MonteCarlo<T>::init(unsigned int const& N_spin, unsigned int const& N_m, Matrice<double> const& H, Array2D<unsigned int> const& sts, Matrice<T> EVec, unsigned int const& nthreads){
	S = new System<T>[nthreads];
	sampling = new std::vector<double>[nthreads];
	E = new double[nthreads];
	err = new double[nthreads];
	for(unsigned int i(0);i<nthreads;i++){
		S[i].init(N_spin,N_m,H,sts,EVec,i);
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

template<typename T>
void MonteCarlo<T>::test_convergence(unsigned int const& thread){
	if( fabs( E[thread] / sampling[thread].size() )  > 1e5 ){ 
		std::cerr<<filename<< " : initial condition lead to a wrong value, restarting the simulation (E="<<E[thread]<<")"<<std::endl;
		E[thread] = 0;
		sampling[thread].clear();
	}  else {
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
		if(cond/m<1e-2 && keep_measuring){ keep_measuring = false; }
	
		output<<S[thread].N_spin
			<<" "<<S[thread].N_m
			<<" "<<sampling[thread].size()
			<<" "<<E[thread] / (sampling[thread].size() * S[thread].N_site)
			<<" "<<err[thread];
		for(unsigned int i(0);i<d.size();i++){
			output<<" "<<d[i];
		}
		output<<Write::endl;
		if(stop.time_limit_reached(time_limit)){
			keep_measuring = false;
			result<<"% the simulation was stopped because it reached the time limit"<<Write::endl;
			std::cerr<<"the simulation was stopped because it reached the time limit"<<std::endl;
		}
	}
}

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

template<typename T>
void MonteCarlo<T>::save(unsigned int const& nthreads){
	for(unsigned int thread(0); thread<nthreads; thread++){
		result<<S[thread].N_spin
			<<" "<<S[thread].N_m
			<<" "<<sampling[thread].size()
			<<" "<<E[thread] / (sampling[thread].size() * S[thread].N_site)
			<<" "<<err[thread]
			<<Write::endl;
	}
}
#endif
