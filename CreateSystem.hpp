#ifndef DEF_CREATESYSTEM
#define DEF_CREATESYSTEM

#include "Parseur.hpp"
#include "Lapack.hpp"
#include "Gnuplot.hpp"

/*!Class that creates a file containing all the necessary information to run a
 * Monte-Carlo simulation.
*/
template<typename Type>
class CreateSystem{
	public:
		/*!Parseur needs N and m, z is the coordination number*/
		CreateSystem(Parseur& P, unsigned int z, std::string filename); 
		/*Simple destructor*/
		virtual ~CreateSystem();

	protected:
		std::string wf_;			//!< type of wavefunction
		unsigned int const m_;		//!< number of unit cell
		unsigned int const N_;		//!< N of SU(N)
		unsigned int const n_;		//!< number of sites
		unsigned int const z_;		//!< coordination number
		double bc_;					//!< boundary condition
		Matrix<unsigned int> sts_;	//!< list of connected sites
		Matrix<Type> T_;			//!< Gutzwiller Hamiltonian
		Matrix<Type> EVec_;			//!< eigenvectors Matrix (transfer Matrix)
		bool successful_;			//!< no degeneracy at the fermi level
		std::string filename_;		//!< filename of the System

		/*!return the neighbours of site i*/
		virtual Vector<unsigned int> get_neighbourg(unsigned int i) = 0;
		/*!compute the array of pairs of swapping sites*/
		void compute_sts();

		/*!compute the eigenvectors from the mean field Hamiltonian*/
		void diagonalize_T(char mat_type);
};

template<typename Type>
CreateSystem<Type>::CreateSystem(Parseur& P, unsigned int z, std::string filename): 
	wf_(P.get<std::string>("wf")),
	m_(P.get<unsigned int>("m")),
	N_(P.get<unsigned int>("N")), 
	n_(N_*m_),
	z_(z),
	bc_(0),
	sts_(n_*z/2,2),
	T_(n_,n_,0.0),
	EVec_(N_*n_,m_),
	successful_(false),
	filename_(filename)
{ }

template<typename Type>
CreateSystem<Type>::~CreateSystem(){}

template<typename Type>
void CreateSystem<Type>::compute_sts(){
	unsigned int k(0);
	Vector<unsigned int> neighbourg;
	for(unsigned int i(0); i<n_;i++){
		neighbourg = get_neighbourg(i);
		for(unsigned int j(0); j<z_/2; j++){
			sts_(k,0) = i;
			sts_(k,1) = neighbourg(j);
			k++;
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
#endif
