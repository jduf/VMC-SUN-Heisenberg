#ifndef DEF_CHAIN
#define DEF_CHAIN

#include "GenericSystem.hpp"

template<typename Type>
class Chain: public GenericSystem<Type>{
	public:
		Chain(unsigned int const& Lx, std::string const& filename);
		virtual ~Chain(){}

	protected:
		Matrix<int> get_neighbourg(unsigned int i) const;
		unsigned int Lx_;//!< number of unit cell along the x-axis
};

template<typename Type>
Chain<Type>::Chain(unsigned int const& Lx, std::string const& filename):
	GenericSystem<Type>(2,filename),
	Lx_(Lx)
{
	std::cout<<"chain"<<std::endl;
	if(Lx_ * this->N_ == this->m_ * this-> n_){ this->compute_links(); }
	else { std::cerr<<"Chain(Lx,filename) : unit cell not compatible with number of site"<<std::endl; }
}

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
