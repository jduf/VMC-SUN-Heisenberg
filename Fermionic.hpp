#ifndef DEF_FERMIONIC
#define DEF_FERMIONIC

#include "Lapack.hpp"
#include "System.hpp"

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
		bool is_degenerate() const { return degenerate_; }

	protected:
		/*!Default Constructor*/
		Fermionic();

		Matrix<Type> T_;	//!< trial Hamiltonian
		Vector<double> eval_;
		Matrix<Type>* EVec_;//!< eigenvectors matrix (transfer matrix)
		bool degenerate_;
		
		/*!compute the eigenvectors from the mean field Hamiltonian*/
		void diagonalize_T();
		void init_fermionic();
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
Fermionic<Type>::Fermionic(Fermionic<Type> const& f):
	System(f),
	T_(f.T_),
	EVec_(f.EVec_?new Matrix<Type>[f.N_]:NULL),
	degenerate_(false)
{
	for(unsigned int i(0);i<N_;i++){ EVec_[i] = f.EVec_[i]; }
}

template<typename Type>
Fermionic<Type>::Fermionic():
	EVec_(NULL),
	degenerate_(false)
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

template<>
inline void Fermionic<double>::diagonalize_T(){
	Lapack<double>(T_,false,'S').eigensystem(eval_,true);
	for(unsigned int c(0);c<N_;c++){
		if(std::abs(eval_(M_(c)) - eval_(M_(c)-1))<1e-12){
			std::cerr<<"Degenerate for the color : "<<c<<std::endl;
			degenerate_= true;
			c=N_;
		}
	}
}

template<>
inline void Fermionic<std::complex<double> >::diagonalize_T(){
	Lapack<std::complex<double> >(T_,false,'H').eigensystem(eval_,true);
	for(unsigned int c(0);c<N_;c++){
		if(are_equal(eval_(M_(c)),eval_(M_(c)-1),1e-12)){
			std::cerr<<"Degenerate for the color : "<<c<<std::endl;
			degenerate_= true;
			c=N_;
		}
	}
}
#endif

