#ifndef DEF_CHAIN
#define DEF_CHAIN

#include "GenericSystem.hpp"

template<typename Type>
class Chain: public GenericSystem<Type>{
	public:
		Chain(std::string const& filename);
		virtual ~Chain();

	protected:
		Matrix<int> get_neighbourg(unsigned int i) const;
		unsigned int a_;//vector of the unit cell
};

template<typename Type>
Chain<Type>::Chain(std::string const& filename):
	GenericSystem<Type>(2,filename),
	a_(this->N_/this->m_)
{
	std::cout<<"chain"<<std::endl;
	unsigned int ppmc(1);
	while(a_*this->m_ != ppmc*this->N_ && a_ < this->n_){ a_ = ++ppmc * this->N_/this->m_; }
	if(a_>this->n_){ std::cerr<<"Chain::Chain(N,n,m,bc,filename) : no unit vector found"<<std::endl; } 
	else { this->compute_links(); }
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
