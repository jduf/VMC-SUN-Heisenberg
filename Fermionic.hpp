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
		/*!Constructor that reads from file*/
		Fermionic(IOFiles& r);
		/*!Destructor*/
		virtual ~Fermionic();
		/*{Forbidden*/
		Fermionic(Fermionic<Type>&&) = delete;
		Fermionic<Type>& operator=(Fermionic<Type>) = delete;
		/*}*/

		Matrix<Type> const& get_EVec() const { return EVec_; }
		Fermionic<Type> const& get_fermionic() const { return (*this); }

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
Fermionic<Type>::Fermionic(Fermionic<Type> const& f):
	System(f),
	same_wf_(f.same_wf_),
	EVec_(f.EVec_?new Matrix<Type>[f.N_]:NULL)
{
	if(same_wf_){ for(unsigned int c(0);c<N_;c++){ EVec_[c] = f.EVec_[0]; } }
	else { for(unsigned int c(0);c<N_;c++){ EVec_[c] = f.EVec_[c]; } }
}

template<typename Type>
Fermionic<Type>::Fermionic(IOFiles& r):
	System(r),
	same_wf_(r.read<bool>()),
	EVec_(N_?new Matrix<Type>[N_]:NULL)
{
	r>>EVec_[0];
	if(same_wf_){ for(unsigned int c(1);c<N_;c++){ EVec_[c] = EVec_[0]; } }
	else { for(unsigned int c(1);c<N_;c++){ r>>EVec_[c]; } }
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

