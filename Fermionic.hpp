#ifndef DEF_FERMIONIC
#define DEF_FERMIONIC

#include "Matrix.hpp"
#include "Vector.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class Fermionic{
	public:
		Fermionic();
		virtual ~Fermionic();

		Matrix<Type> const& get_EVec() const { return EVec_;}

	protected:
		Matrix<Type> T_;			//!< Gutzwiller Hamiltonian
		Matrix<Type> EVec_;			//!< eigenvectors Matrix (transfer Matrix)
		
	private:
		Fermionic(Fermionic<Type> const& f);
};

/*constructors and destructor and initialization*/
/*{*/
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

