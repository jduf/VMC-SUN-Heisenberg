#ifndef DEF_FERMIONIC
#define DEF_FERMIONIC

#include "System.hpp"
#include "Lapack.hpp"

/*{*//*!Class that connects GenericSystem and MCSystem

	   This class is essential because it allows to transfer relevant variables
	   from an instance of a child of GenericSystem (e.g. ChainPolymerized) to
	   an instance of a child of MCSystem (e.g. SystemFermionic)

	   The most important variable is EVec_ which contains, for a given color
	   c, the M_(c) lowest eigenvectors of the trial Hamiltonian.*//*}*/
template<typename Type>
class Fermionic : public virtual System{
	public:
		/*!Copy constructor*/
		Fermionic(Fermionic<Type> const& F);
		/*!Constructor that reads from file*/
		Fermionic(IOFiles& r);
		/*!Destructor*/
		virtual ~Fermionic();
		/*{Forbidden*/
		Fermionic(Fermionic<Type>&&) = delete;
		Fermionic<Type>& operator=(Fermionic<Type>) = delete;
		/*}*/

		Fermionic<Type> const& get_Fermionic() const { return (*this); }

		Matrix<Type> const& get_EVec(unsigned int const& c) const { assert(c<N_); return EVec_[c]; }
		Matrix<Type>* get_EVec() const { return EVec_; }
		bool const& use_same_wf() const { return same_wf_; }

	protected:
		/*!Constructor*/
		Fermionic();

		bool same_wf_;		//!< true if the same wavefunction is used for all colors
		Matrix<Type>* EVec_;//!< eigenvectors matrix (transfer matrix)

		/*!Sets the eigenvectors' matrices*/
		void init_fermionic();
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
Fermionic<Type>::Fermionic(Fermionic<Type> const& F):
	System(F),
	same_wf_(F.same_wf_),
	EVec_(new Matrix<Type>[N_])
{
	if(F.EVec_){
		if(same_wf_){ for(unsigned int c(0);c<N_;c++){ EVec_[c] = F.EVec_[0]; } }
		else { for(unsigned int c(0);c<N_;c++){ EVec_[c] = F.EVec_[c]; } }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : EVec_ == NULL"<<std::endl; }
}

template<typename Type>
Fermionic<Type>::Fermionic(IOFiles& r):
	System(r),
	same_wf_(r.read<bool>()),
	EVec_(new Matrix<Type>[N_])
{
	if(N_){
		r>>EVec_[0];
		if(same_wf_){ for(unsigned int c(1);c<N_;c++){ EVec_[c] = EVec_[0]; } }
		else { for(unsigned int c(1);c<N_;c++){ r>>EVec_[c]; } }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : N_ == 0"<<std::endl; }
}

template<typename Type>
Fermionic<Type>::Fermionic():
	same_wf_(true),
	EVec_(NULL)
{}

template<typename Type>
void Fermionic<Type>::init_fermionic(){
	if(!EVec_){ EVec_ = new Matrix<Type>[N_]; }
	for(unsigned int c(0);c<N_;c++){ EVec_[c].set(n_,M_(c)); }
}

template<typename Type>
Fermionic<Type>::~Fermionic(){
	if(EVec_){ delete[] EVec_; }
}
/*}*/
#endif
