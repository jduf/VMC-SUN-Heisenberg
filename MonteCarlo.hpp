#ifndef DEF_MONTECARLO
#define DEF_MONTECARLO

#include "System.hpp"
#include "Array2D.hpp"
#include "Write.hpp"

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
		void test(unsigned int const& thread);
		void init(unsigned int N_spin, unsigned int N_m, Matrice<double> const& H, Array2D<unsigned int> const& sts, Matrice<T> EVec, unsigned int nthreads);

	private:
		/*!Forbids the copy constructor*/
		MonteCarlo(MonteCarlo const& mc);
		/*!Forbids the assignment operator*/
		MonteCarlo const& operator=(MonteCarlo const& mc);

		bool binning_analysis(unsigned int& iter, unsigned int const& thread);

		/*!Computes the mean of a std::vector<double> */
		double mean(std::vector<double> const& v);

		/*!Compute the variance of a std::vector<double> */
		double delta(std::vector<double> const& v, double m);
	
		System<T>* S; //!< Pointer to a system 
		std::vector<double>* sampling; //<! Stores all the values that MC considers
		double* E; //!< Value that the MC algorithm tries to compute
		unsigned int const N_MC; //!< Number of measure to do before doing a binning analyis
		Write output; //!< Text file of name "filename" that stores the output of the MC algorithm
		std::string filename; //!< Name of the output file
};

template<typename T>
MonteCarlo<T>::MonteCarlo(std::string filename):
	N_MC(1e4),
	output(filename+".out"),
	filename(filename)
{}

template<typename T>
MonteCarlo<T>::~MonteCarlo(){}

template<typename T>
void MonteCarlo<T>::init(unsigned int N_spin, unsigned int N_m, Matrice<double> const& H, Array2D<unsigned int> const& sts, Matrice<T> EVec, unsigned int nthreads){
	std::cout<<nthreads<<std::endl;
	S = new System<T>[nthreads];
	sampling = new std::vector<double>[nthreads];
	E = new double[nthreads];
	for(unsigned int i(0);i<nthreads;i++){
		S[i].init(N_spin,N_m,H,sts,EVec,i);
		//S[i].swap();
		//S[i].update();
		//S[i].print();
	}
}

template<typename T>
bool MonteCarlo<T>::binning_analysis(unsigned int& iter, unsigned int const& thread){
	if(iter==N_MC){
		std::vector<double> bin(sampling[thread]);
		std::vector<double> d(0);
		unsigned int l(0);
		while(pow(2,l+1)*30<sampling[thread].size()){
			l++;
			std::vector<double> bin2(bin.size()/2);
			for(unsigned int i(0);i<bin2.size();i++){
				bin2[i] = (bin[2*i] + bin[2*i+1])/2.0;
			}
			d.push_back(delta(bin,mean(bin)));
			bin.resize(bin2.size());
			bin=bin2;
		}

		output<<S[thread].N_spin
		 <<" "<<S[thread].N_m
		 <<" "<<sampling[thread].size()
		 <<" "<<E[thread] / (sampling[thread].size() * S[thread].N_site);
		for(unsigned int i(0);i<d.size();i++){
			output<<" "<<d[i];
		}
		output<<Write::endl;
		if( std::abs( E[thread] / (sampling[thread].size() * S[thread].N_site) )  > 100 ){ 
			std::cerr<<filename<< " : initial condition lead to a wrong value, resetting the simulation"<<std::endl;
			E[thread] = 0;
			sampling[thread].clear();
		} 
		if(fabs(d[d.size()-1]) > 1e-2){iter=0; return true;}
		else {
			output<<S[thread].N_spin
				<<" "<<S[thread].N_m
				<<" "<<sampling[thread].size()
				<<" "<<E[thread] / (sampling[thread].size() * S[thread].N_site)
				<<" "<<d[d.size()-1]
				<<Write::endl;
			return false; 
		}
	} else {
		return true;
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
double MonteCarlo<T>::delta(std::vector<double> const& v,double m){
	double d(0.0);
	unsigned int N(v.size());
	for(unsigned int i(0);i<N;i++){
		d += (v[i]-m)*(v[i]-m);
	}
	return sqrt(d / (N*(N-1))); 
}

template<typename T>
void MonteCarlo<T>::test(unsigned int const& thread){
	S[thread].print();
	S[thread].swap();
	std::cout<<"updating : ratio="<<S[thread].ratio()<<std::endl;
	S[thread].update();
	std::cout<<S[thread].compute_energy()<<std::endl;
	S[thread].print();
}
#endif
