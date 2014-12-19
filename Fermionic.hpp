#ifndef DEF_FERMIONIC
#define DEF_FERMIONIC

#include "System.hpp"
#include "Lapack.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class Fermionic : public virtual System{
	public:
		/*!Copy Constructor*/
		Fermionic(Fermionic<Type> const& f);
		/*!Destructor*/
		virtual ~Fermionic();

		Matrix<Type> const& get_EVec() const { return EVec_; }
		Fermionic<Type> const& get_fermionic() const { return (*this); }

	protected:
		/*!Default Constructor*/
		Fermionic();

		Matrix<Type>* EVec_;//!< eigenvectors matrix (transfer matrix)
		
		/*!compute the eigenvectors from the mean field Hamiltonian*/
		void init_fermionic();
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
Fermionic<Type>::Fermionic(Fermionic<Type> const& f):
	System(f),
	EVec_(f.EVec_?new Matrix<Type>[f.N_]:NULL)
{
	for(unsigned int c(0);c<N_;c++){ EVec_[c] = f.EVec_[c]; }
}

template<typename Type>
Fermionic<Type>::Fermionic():
	EVec_(NULL)
{}

template<typename Type>
void Fermionic<Type>::init_fermionic(){
	if(!EVec_){ EVec_ = new Matrix<Type>[N_];}
	for(unsigned int c(0);c<N_;c++){ EVec_[c].set(n_,M_(c)); }
}

template<typename Type>
Fermionic<Type>::~Fermionic(){
	if(EVec_){ delete[] EVec_; }
}
/*}*/
#endif

