#ifndef DEF_MONTECARLO
#define DEF_MONTECARLO

#include "System.hpp"
#include "Array2D.hpp"
#include "Write.hpp"

#include <vector>

template <typename T>
class MonteCarlo{
	public:
		MonteCarlo(System<T>* S, Matrice<double> const& H, Array2D<unsigned int> const& sts, std::string filename);
		~MonteCarlo();

		void run();

	private:
		MonteCarlo(MonteCarlo const& MonteCarlo);
		MonteCarlo const& operator=(MonteCarlo const& mc);

		System<T>* S;
		Matrice<double> const H;
		Array2D<unsigned int> const sts;
		std::vector<double> sampling;
		double E;
		unsigned int const N_MC;
		Write output;
		std::string filename;

		void compute_energy();
		bool binning_analyse(unsigned int& iter);
		double mean(std::vector<double> v);
		double delta(std::vector<double> v, double m);
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
bool MonteCarlo<T>::binning_analyse(unsigned int& iter){
	if(iter==N_MC){
		std::vector<double> bin(sampling);
		std::vector<double> d(0);
		unsigned int l(0);
		while(pow(2,l+1)*30<N_MC){
			l++;
			std::vector<double> bin2(bin.size()/2);
			for(unsigned int i(0);i<bin2.size();i++){
				bin2[i] = (bin[2*i] + bin[2*i+1])/2.0;
			}
			d.push_back(delta(bin,mean(bin)));
			bin.resize(bin2.size());
			bin=bin2;
		}

		output<<" "<<S->N_spin
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
		if(fabs(d[d.size()-1]) > 5e-2){iter=0; return true;}
		else { return false; }
	} else {
		return true;
	}
}

template<typename T>
double MonteCarlo<T>::mean(std::vector<double> v){
	double m(0.0);
	for(unsigned int i(0);i<v.size();i++){
		m += v[i];
	}
	return m/v.size();
}

template<typename T>
double MonteCarlo<T>::delta(std::vector<double> v,double m){
	double d(0.0);
	for(unsigned int i(0);i<v.size();i++){
		d += (v[i]-m)*(v[i]-m);
	}
	return sqrt(d / (v.size()*(v.size()-1))); 
}
#endif
