#ifndef DEF_FERMIONIC
#define DEF_FERMIONIC

#include "Matrix.hpp"
#include "Vector.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class Fermionic{
	public:
		Fermionic(Fermionic<Type> const& f);
		Fermionic();
		virtual ~Fermionic();

		Matrix<Type> const& get_EVec() const { return EVec_;}

		Fermionic<Type> const& get_fermionic() const { return (*this);}

	protected:
		Matrix<Type> T_;			//!< Gutzwiller Hamiltonian
		Matrix<Type> EVec_;			//!< eigenvectors Matrix (transfer Matrix)
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
Fermionic<Type>::Fermionic(Fermionic<Type> const& f):
	T_(f.T_),
	EVec_(f.EVec_)
{
	std::cout<<"ok normal Fermionic"<<std::endl;
}

template<typename Type>
Fermionic<Type>::Fermionic() {
	std::cout<<"ok default Fermionic"<<std::endl;
}

template<typename Type>
Fermionic<Type>::~Fermionic(){
	std::cout<<"~Fermionic"<<std::endl;
}
/*}*/
#endif

