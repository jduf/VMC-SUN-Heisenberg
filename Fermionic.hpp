#ifndef DEF_FERMIONIC
#define DEF_FERMIONIC

#include "Lapack.hpp"
#include "System.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class Fermionic : public virtual System{
	public:
		Fermionic(Fermionic<Type> const& f);
		virtual ~Fermionic();

		Matrix<Type> const& get_EVec() const { return EVec_;}

		Fermionic<Type> const& get_fermionic() const { return (*this);}
		bool is_degenerate() const { return degenerate_;}

	protected:
		Matrix<Type> T_;	//!< trial Hamiltonian
		Vector<double> eval_;
		Matrix<Type>* EVec_;//!< eigenvectors matrix (transfer matrix)
		bool degenerate_;
		
		/*!compute the eigenvectors from the mean field Hamiltonian*/
		void diagonalize_T();
		bool is_degenerate(unsigned int const& c);
		void init_fermionic();

		Fermionic();
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
Fermionic<Type>::Fermionic(Fermionic<Type> const& f):
	System(f),
	T_(f.T_),
	EVec_(NULL),
	degenerate_(false)
{
	std::cerr<<"Fermionic copy do something for the partial copy"<<std::endl;
	if(f.EVec_){
		EVec_ = new Matrix<Type>[f.N_];
		for(unsigned int i(0);i<N_;i++){ EVec_[i] = f.EVec_[i]; }
	}
}

template<typename Type>
Fermionic<Type>::Fermionic():
	EVec_(NULL),
	degenerate_(false)
{std::cout<<"fermionic default"<<N_<<std::endl;}

template<typename Type>
Fermionic<Type>::~Fermionic(){
	std::cout<<"destroy Fermionic"<<std::endl;
	if(EVec_){ delete[] EVec_; }
}
/*}*/

template<>
inline void Fermionic<double>::diagonalize_T(){
	Lapack<double> ES(&T_,false,'S');
	ES.eigensystem(&eval_,true);
}

template<>
inline void Fermionic<std::complex<double> >::diagonalize_T(){
	Lapack<std::complex<double> >ES(&T_,false,'H');
	ES.eigensystem(&eval_,true);
}

template<typename Type>
bool Fermionic<Type>::is_degenerate(unsigned int const& c) {
	if(std::abs(eval_(M_(c)) - eval_(M_(c)-1))<1e-12){
		std::cerr<<"Fermi level degenerate"<<std::endl;
		//degenerate_=true;
	}
	return false;
		//return true;
		//} else {
		//return false;
}

template<typename Type>
void Fermionic<Type>::init_fermionic(){
	if(!EVec_){ EVec_ = new Matrix<Type>[N_];}
	for(unsigned int c(0);c<N_;c++){ EVec_[c].set(n_,M_(c)); }
}
#endif

