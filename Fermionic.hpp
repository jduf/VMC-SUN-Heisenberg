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

		Matrix<Type>* EVec_;//!< eigenvectors matrix (transfer matrix)
		bool degenerate_;
		
		/*!compute the eigenvectors from the mean field Hamiltonian*/
		void init_fermionic();
		void diagonalize_H(Matrix<Type>& H);
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
Fermionic<Type>::Fermionic(Fermionic<Type> const& f):
	System(f),
	EVec_(f.EVec_?new Matrix<Type>[f.N_]:NULL),
	degenerate_(f.degenerate_)
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
inline void Fermionic<double>::diagonalize_H(Matrix<double>& H){
	Vector<double> eval;
	Lapack<double>(H,false,'S').eigensystem(eval,true);
	for(unsigned int c(0);c<N_;c++){
		if(std::abs(eval(M_(c)) - eval(M_(c)-1))<1e-12){
			std::cerr<<"Degenerate for the color : "<<c<<std::endl;
			degenerate_= true;
			c=N_;
		}
	}
}

template<>
inline void Fermionic<std::complex<double> >::diagonalize_H(Matrix<std::complex<double> >& H){
	Vector<double> eval;
	Lapack<std::complex<double> >(H,false,'H').eigensystem(eval,true);
	for(unsigned int c(0);c<N_;c++){
		if(are_equal(eval(M_(c)),eval(M_(c)-1),1e-12)){
			std::cerr<<"Degenerate for the color : "<<c<<std::endl;
			degenerate_= true;
			c=N_;
		}
	}
}

#endif

