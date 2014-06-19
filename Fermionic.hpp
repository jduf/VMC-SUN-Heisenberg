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
		bool is_degenerate() const { return degenerate_; }

	protected:
		Matrix<Type> T_;	//!< trial Hamiltonian
		Matrix<Type> EVec_;	//!< eigenvectors matrix (transfer matrix)
		bool degenerate_;	//!< no degeneracy at the fermi level
		
		/*!compute the eigenvectors from the mean field Hamiltonian*/
		void diagonalize_T(char mat_type);

		Fermionic(){std::cout<<"fermionic default"<<std::endl;};
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
Fermionic<Type>::Fermionic(Fermionic<Type> const& f):
	System(f),
	T_(f.T_),
	EVec_(f.EVec_),
	degenerate_(false)
{std::cout<<"fermionic copy"<<std::endl;}

template<typename Type>
Fermionic<Type>::~Fermionic(){}
/*}*/

template<typename Type>
void Fermionic<Type>::diagonalize_T(char mat_type){
	std::cout<<"ok"<<T_<<std::endl;
	Lapack<Type> ES(&T_,false, mat_type);
	Vector<double> EVal;
	ES.eigensystem(&EVal,true);
	if(std::abs(EVal(M_) - EVal(M_-1))<1e-12){
		std::cerr<<"Fermi level degenerate"<<std::endl;
		degenerate_ = true;
	} else {
		degenerate_ = false;
	}
}
#endif

