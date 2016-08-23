#ifndef DEF_BOSONIC
#define DEF_BOSONIC

#include "System.hpp"

/*!Class that contains the information on the state*/
template<typename Type>
class Bosonic : public virtual System{
	public:
		/*!Copy Constructor*/
		Bosonic(Bosonic<Type> const& b);
		/*!Constructor that reads from file*/
		Bosonic(IOFiles& r);
		/*!Default destructor*/
		virtual ~Bosonic() = default;
		/*{Forbidden*/
		Bosonic(Bosonic<Type>&&) = delete;
		Bosonic<Type>& operator=(Bosonic<Type>) = delete;
		/*}*/

		Vector<unsigned int> const& get_sl() const { return sl_; }
		Matrix<unsigned int> const& get_nn() const { return nn_; }
		Matrix<unsigned int> const& get_cc() const { return cc_; }
		Matrix<double> const& get_nu() const { return nu_; }
		Matrix<Type> const& get_omega() const { return omega_; }

	protected:
		/*!Default Constructor*/
		Bosonic() = default;

		Vector<unsigned int> sl_;//!< sl_(i) is the sublattice to which the site i belongs to
		Matrix<unsigned int> nn_;//!< nn_(i,j):jth neighbour of the ith site
		Matrix<unsigned int> cc_;//!< connect a combination of two color to one nu_
		Matrix<double> nu_;      //!< nu_(i,j): i link and j factor
		Matrix<Type> omega_;     //!< \warning maybe not a complex<double> but a double. before it was Type

		void init_bosonic(unsigned int const& z, Matrix<double> const& nu);
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
Bosonic<Type>::Bosonic(IOFiles& r):
	System(r),
	sl_(r),
	nn_(r),
	cc_(r),
	nu_(r),
	omega_(r)
{}

template<typename Type>
void Bosonic<Type>::init_bosonic(unsigned int const& z, Matrix<double> const& nu){
	nu_ = nu;
	nn_.set(n_,z);
	cc_.set(N_,N_);
	sl_.set(n_);
	omega_.set(N_,N_,1.0);
}
/*}*/
#endif
