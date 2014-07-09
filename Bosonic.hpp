#ifndef DEF_BOSONIC
#define DEF_BOSONIC

#include "System.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class Bosonic : public virtual System{
	public:
		Bosonic(Bosonic<Type> const& b);
		virtual ~Bosonic(){}

		Vector<unsigned int> const& get_sl() const { return sl_;}
		Matrix<unsigned int> const& get_nn() const { return nn_;}
		Matrix<unsigned int> const& get_cc() const { return cc_;}
		Matrix<double> const& get_nu() const { return nu_;}
		Matrix<Type> const& get_omega() const { return omega_;}

		Bosonic<Type> const& get_bosonic(){ return (*this);}

	protected:
		Bosonic(){}

		Vector<unsigned int> sl_;
		Matrix<unsigned int> nn_; //!< nn_(i,j):jth neighbour of the ith site
		Matrix<unsigned int> cc_;
		Matrix<double> nu_;
		Matrix<Type> omega_; //!< \warning maybe not a complex<double> but a double. before it was Type

		void init_bosonic(unsigned int const& z);
};

/*constructors and destructor*/
/*{*/
template<typename Type>
Bosonic<Type>::Bosonic(Bosonic<Type> const& b):
	System(b),
	sl_(b.sl_),
	nn_(b.nn_),
	cc_(b.cc_),
	nu_(b.nu_),
	omega_(b.omega_)
{}

template<typename Type>
void Bosonic<Type>::init_bosonic(unsigned int const& z){
	nn_.set(n_,z);
	cc_.set(N_,N_);
	sl_.set(n_);
	omega_.set(N_,N_,1.0);
}
/*}*/
#endif

