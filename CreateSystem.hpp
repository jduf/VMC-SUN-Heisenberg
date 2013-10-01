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
		/*!Parseur needs N and m, z is the coordination number*/
		CreateSystem(Parseur& P, unsigned int z); 
		~CreateSystem();

	protected:
		std::string wf_;			//!< type of wavefunction
		unsigned int const m_;		//!< number of unit cell
		unsigned int const N_;		//!< N of SU(N)
		unsigned int const n_;		//!< number of sites
		double bc_;					//!< boundary condition
		Matrix<unsigned int> sts_;	//!< list of connected sites
		Matrix<int> H_;				//!< SU(N) Matrix
		Matrix<Type> T_;			//!< Gutzwiller Hamiltonian
		Matrix<Type> EVec_;			//!< eigenvectors Matrix (transfer Matrix)
		bool successful_;			//!< no degeneracy at the fermi level
		bool study_system_;		//!< no degeneracy at the fermi level

		/*!compute the eigenvectors from the mean field Hamiltonian*/
		void diagonalize_T(char mat_type);
		/*!compute the array of pairs of swapping sites*/
		void compute_sts();
		/*!Evaluate the value of an operator O as <bra|O|ket>*/
		std::complex<double> projection(Matrix<double> const& O, Matrix<std::complex<double> > const& base, unsigned int bra, unsigned int ket);
};

template<typename Type>
CreateSystem<Type>::CreateSystem(Parseur& P, unsigned int z): 
	wf_(P.get<std::string>("wf")),
	m_(P.get<unsigned int>("m")),
	N_(P.get<unsigned int>("N")), 
	n_(N_*m_),
	bc_(0),
	sts_(n_*z/2,2),
	H_(n_,n_,0.0),
	T_(n_,n_,0.0),
	EVec_(N_*n_,m_),
	successful_(false),
	study_system_(false)
{
	P.set("study",study_system_);
}

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
void CreateSystem<Type>::diagonalize_T(char mat_type){
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
