#ifndef DEF_CREATESYSTEM
#define DEF_CREATESYSTEM

#include "Parseur.hpp"
#include "Lapack.hpp"
#include "Gnuplot.hpp"

/*!Class that creates a file containing all the necessary information to run a
 * Monte-Carlo simulation.
 *
 *  
 *
*/
template<typename Type>
class CreateSystem{
	public:
		CreateSystem(Parseur& P,unsigned int N_n);
		~CreateSystem();

	protected:
		unsigned int const m_, N_, n_;//!<
		double bc_;				//!<
		Matrix<unsigned int> sts_;//!< 
		Matrix<int> H_;			//!< H_opping Matrix
		Matrix<Type> T_;			//!< GuT_zwiller H_amilT_onian
		Matrix<Type> EVec_;		//!< eigenVectors Matrix (transfer Matrix)
		bool successful_;		//!< no degeneracy aT_ T_H_e fermi level

		/*!compute T_H_e eigenVectors from T_H_e mean field H_amilT_onian*/
		void diagonalize_EVec(char mat_type);
		/*!compute T_H_e array of pairs of swapping siT_es*/
		void compute_sts();

		std::complex<double> projection(Matrix<double> const& O, Matrix<std::complex<double> > const& base, unsigned int bra, unsigned int ket);
};

template<typename Type>
CreateSystem<Type>::CreateSystem(Parseur& P, unsigned int N_n):
	m_(P.get<unsigned int>("N_m")),
	N_(P.get<unsigned int>("N_spin")), 
	n_(N_*m_),
	bc_(0),
	sts_(n_*N_n/2,2),
	H_(n_,n_,0.0),
	T_(n_,n_,0.0),
	EVec_(N_*n_,m_),
	successful_(false)
{ }

template<typename Type>
CreateSystem<Type>::~CreateSystem(){ }

template<typename Type>
void CreateSystem<Type>::compute_sts(){
	unsigned int k(0);
	for(unsigned int i(0); i<n_;i++){
		for(unsigned int j(i+1); j<n_;j++){
			if ( H_(i,j) != 0){
				sts_(k,0) = i;
				sts_(k,1) = j;
				k++;
			}
		}
	}
}

template<typename Type>
void CreateSystem<Type>::diagonalize_EVec(char mat_type){
	Lapack<Type> ES(&T_,false, mat_type);
	Vector<double> EVal;
	ES.eigensystem(&EVal,true);
	//std::cout<<EVal<<std::endl;
	if(std::abs(EVal(m_) - EVal(m_-1))>1e-10){ successful_ = true; }
}

template<typename Type>
std::complex<double> CreateSystem<Type>::projection(Matrix<double> const& O, Matrix<std::complex<double> > const& base, unsigned int bra, unsigned int ket){
	Vector<std::complex<double> > tmp(n_,0.0);
	std::complex<double> out(0.0);
	for(unsigned int i(0);i<n_;i++){
		for(unsigned int j(0);j<n_;j++){
			tmp(i) += O(i,j)*base(j,ket);
		}
	}
	for(unsigned int i(0);i<n_;i++){
		out += tmp(i)*std::conj(base(i,bra));
	}
	return out;
}
#endif
