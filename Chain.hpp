#ifndef DEF_CHAIN
#define DEF_CHAIN

#include "GenericSystem.hpp"

template<typename Type>
class Chain: public GenericSystem<Type>{
	public:
		Chain(unsigned int N, unsigned int n, unsigned int m, int bc, std::string filename);
		virtual ~Chain();

	protected:
		Matrix<int> get_neighbourg(unsigned int i) const;
		unsigned int a_;//vector of the unit cell
};

template<typename Type>
Chain<Type>::Chain(unsigned int N, unsigned int n, unsigned int m, int bc, std::string filename):
	GenericSystem<Type>(N,n,m,bc,2,filename),
	a_(N/m)
{
	if(a_*m != N){ std::cerr<<"Chain::Chain(N,n,m,z,filename) : The unit cell is not well defined"<<std::endl; }
	this->compute_sts();
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
