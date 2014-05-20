#ifndef DEF_FERMIONIC
#define DEF_FERMIONIC

#include "Matrix.hpp"
#include "Vector.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class Fermionic{
	public:
		Fermionic();
		Fermionic(Fermionic const& f);
		virtual ~Fermionic();

		Matrix<Type> get_EVec() const { return EVec_;}
		Fermionic<Type>* get_fermionic(){ return this;}

	protected:
		Matrix<Type> T_;			//!< Gutzwiller Hamiltonian
		Matrix<Type> EVec_;			//!< eigenvectors Matrix (transfer Matrix)
};

/*constructors and destructor and initialization*/
/*{*/
template<typename Type>
Fermionic<Type>::Fermionic()
{
	std::cout<<"ok Fermionic"<<std::endl;
}

template<typename Type>
Fermionic<Type>::Fermionic(Fermionic const& f):
	T_(f.T_),
	EVec_(f.EVec_)
{ 
	std::cout<<"ok copy Fermionic"<<std::endl;
}

template<typename Type>
Fermionic<Type>::~Fermionic(){}
/*}*/
#endif

