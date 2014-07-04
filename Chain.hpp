#ifndef DEF_CHAIN
#define DEF_CHAIN

#include "GenericSystem.hpp"

template<typename Type>
class Chain: public GenericSystem<Type>{
	public:
		Chain(unsigned int const& spuc, std::string const& filename);
		virtual ~Chain()=0;

	protected:
		Matrix<int> get_neighbourg(unsigned int i) const;
		unsigned int Lx_;//!< number of unit cell along the x-axis
		unsigned int spuc_;//!< site per unit cell
};

template<typename Type>
Chain<Type>::Chain(unsigned int const& spuc, std::string const& filename):
	GenericSystem<Type>(2,filename),
	Lx_(this->n_/spuc),
	spuc_(spuc)
{
	this->compute_links(); 
	this->status_--;
}

template<typename Type>
Chain<Type>::~Chain(){}

template<typename Type>
Matrix<int> Chain<Type>::get_neighbourg(unsigned int i) const {
	Matrix<int> nb(this->z_,2,1);
	if( i != this->n_-1){ nb(0,0) = i+1;}
	else {
		nb(0,0) = 0;
		nb(0,1) = this->bc_;
	}
	if( i != 0){ nb(1,0) = i-1;}
	else {
		nb(1,0) = this->n_-1;
		nb(1,1) = this->bc_;
	}
	return nb;
}
#endif
