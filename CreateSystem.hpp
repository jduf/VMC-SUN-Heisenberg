#ifndef DEF_CREATESYSTEM
#define DEF_CREATESYSTEM

#include "Parseur.hpp"
#include "Write.hpp"
#include "Lapack.hpp"

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
		unsigned int const N_m, N_spin, N_site;//!<
		double bc;				//!<
		Matrix<unsigned int> sts;//!< 
		Matrix<int> H;			//!< hopping matrix
		Matrix<Type> T;			//!< Gutzwiller Hamiltonian
		Matrix<Type> EVec;		//!< eigenvectors matrix (transfer matrix)
		bool successful;		//!< no degeneracy at the fermi level
		char mat_type;			//!< matrix type

		/*!Compute the eigenvectors from the mean field hamiltonian*/
		void compute_EVec();
		/*!Compute the array of pairs of swapping sites*/
		void compute_sts();
};

template<typename Type>
CreateSystem<Type>::CreateSystem(Parseur& P, unsigned int N_n):
	N_m(P.get<unsigned int>("N_m")),
	N_spin(P.get<unsigned int>("N_spin")), 
	N_site(N_spin*N_m),
	bc(0),
	sts(N_spin*N_m*N_n/2,2),
	H(N_spin*N_m,N_spin*N_m,0.0),
	T(N_spin*N_m,N_spin*N_m,0.0),
	EVec(N_spin*N_spin*N_m,N_m),
	successful(false),
	mat_type('U')
{ }

template<typename Type>
CreateSystem<Type>::~CreateSystem(){ }

template<typename Type>
void CreateSystem<Type>::compute_sts(){
	unsigned int k(0);
	for(unsigned int i(0); i<N_site;i++){
		for(unsigned int j(i+1); j<N_site;j++){
			if ( H(i,j) != 0){
				sts(k,0) = i;
				sts(k,1) = j;
				k++;
			}
		}
	}
}

template<typename Type>
void CreateSystem<Type>::compute_EVec(){
	Lapack<Type> ES(&T,false, mat_type);
	Matrix<double> EVal;
	ES.eigensystem(EVal);
	if(std::abs(EVal(N_m) - EVal(N_m-1))>1e-10){ successful = true; }
}
#endif
