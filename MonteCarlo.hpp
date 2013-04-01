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
		//{Description
		/*!Constructors */
		//}
		MonteCarlo(System<T>* S, Matrice<double> const& H, Array2D<unsigned int> const& sts, std::string filename); 
		~MonteCarlo();

		//{Description
		/*!Forbids the copy constructor*/
		//}
		void run();

	private:
		//{Description
		/*!Forbids the copy constructor*/
		//}
		MonteCarlo(MonteCarlo const& MonteCarlo);
		//{Description
		/*!Forbids the assignment operator*/
		//}
		MonteCarlo const& operator=(MonteCarlo const& mc);

		//{Description
		/*!Computes the energy of the configuration
		 *
		 * This function is called after N_site swap() in order to decorrelate
		 * two successive configurations
		 * */
		//}
		void compute_energy();

		bool binning_analysis(unsigned int& iter);

		//{Description
		/*!Computes the mean of a std::vector<double>
		 *
		 * Used during the binning_analysis()
		 * */
		//}
		double mean(std::vector<double> const& v);

		//{Description
		/*!Compute the variance of a std::vector<double>
		 *
		 * Used during the binning_analysis()
		 */
		//}
		double delta(std::vector<double> const& v, double m);
	
		System<T>* S; //!< Pointer to a system 
		Matrice<double> const H; //!< Hamiltonian that is used to compute the energy (might not be the best place to store this value)
		Array2D<unsigned int> const sts;//!< Arrray that is used to compute the energy (might not be the best place to store this value)
		std::vector<double> sampling; //<! Stores all the values that MC considers
		double E; //!< Value that the MC algorithm tries to compute
		unsigned int const N_MC; //!< Number of measure to do before doing a binning analyis
		Write output; //!< Text file of name "filename" that stores the output of the MC algorithm
		std::string filename; //!< Name of the output file

};

template<typename T>
MonteCarlo<T>::MonteCarlo(System<T>* S, Matrice<double> const& H, Array2D<unsigned int> const& sts,std::string filename):
	S(S),
	H(H),
	sts(sts),
	sampling(0),
	E(0.0),
	N_MC(1e4),
	output(filename+".out"),
	filename(filename)
{}

template<typename T>
MonteCarlo<T>::~MonteCarlo(){}

template<typename T>
bool MonteCarlo<T>::binning_analysis(unsigned int& iter){
	if(iter==N_MC){
		std::vector<double> bin(sampling);
		std::vector<double> d(0);
		unsigned int l(0);
		while(pow(2,l+1)*30<sampling.size()){
			l++;
			std::vector<double> bin2(bin.size()/2);
			for(unsigned int i(0);i<bin2.size();i++){
				bin2[i] = (bin[2*i] + bin[2*i+1])/2.0;
			}
			d.push_back(delta(bin,mean(bin)));
			bin.resize(bin2.size());
			bin=bin2;
		}

		output<<S->N_spin
		 <<" "<<S->N_m
		 <<" "<<sampling.size()
		 <<" "<<E / (sampling.size() * S->N_site);
		for(unsigned int i(0);i<d.size();i++){
			output<<" "<<d[i];
		}
		output<<Write::endl;
		if( std::abs( E / (sampling.size() * S->N_site) )  > 100 ){ 
			std::cerr<<filename<< " : initial condition lead to a wrong value, resetting the simulation"<<std::endl;
			E = 0;
			sampling.clear();
		} 
		if(fabs(d[d.size()-1]) > 1e-2){iter=0; return true;}
		else {
			output<<S->N_spin
				<<" "<<S->N_m
				<<" "<<sampling.size()
				<<" "<<E / (sampling.size() * S->N_site)
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
#endif
