#ifndef DEF_BOSONIC
#define DEF_BOSONIC

#include "Matrix.hpp"
#include "Vector.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class Bosonic{
	public:
		Bosonic();
		virtual ~Bosonic();

		Vector<unsigned int> const& get_sl() const { return sl_;}
		Matrix<unsigned int> const& get_nn() const { return nn_;}
		Matrix<unsigned int> const& get_cc() const { return cc_;}
		Matrix<double> const& get_nu() const { return nu_;}
		Matrix<Type> const& get_omega() const { return omega_;}

	protected:
		Vector<unsigned int> sl_;
		Matrix<unsigned int> nn_; //!< nn_(i,j):jth neighbour of the ith site
		Matrix<unsigned int> cc_;
		Matrix<double> nu_;
		Matrix<Type> omega_;

	private:
		Bosonic(Bosonic<Type> const& b);
};

/*constructors and destructor*/
/*{*/
template<typename Type>
Bosonic<Type>::Bosonic(){}

template<typename Type>
Bosonic<Type>::~Bosonic(){}
/*}*/
#endif

