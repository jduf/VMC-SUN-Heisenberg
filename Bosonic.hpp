#ifndef DEF_BOSONIC
#define DEF_BOSONIC

#include "System.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class Bosonic : public virtual System{
	public:
		Bosonic(Bosonic<Type> const& b);
		virtual ~Bosonic();

		Vector<unsigned int> const& get_sl() const { return sl_;}
		Matrix<unsigned int> const& get_nn() const { return nn_;}
		Matrix<unsigned int> const& get_cc() const { return cc_;}
		Matrix<double> const& get_nu() const { return nu_;}
		Matrix<Type> const& get_omega() const { return omega_;}

		Bosonic<Type> const& get_bosonic(){ return (*this);}

	protected:
		Vector<unsigned int> sl_;
		Matrix<unsigned int> nn_; //!< nn_(i,j):jth neighbour of the ith site
		Matrix<unsigned int> cc_;
		Matrix<double> nu_;
		Matrix<Type> omega_;

		Bosonic(){std::cout<<"bosonic default"<<std::endl;};
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
Bosonic<Type>::~Bosonic(){}
/*}*/
#endif

