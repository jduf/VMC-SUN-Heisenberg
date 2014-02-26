#ifndef DEF_CHAIN
#define DEF_CHAIN

#include "GenericSystem.hpp"

template<typename Type>
class Chain: public GenericSystem<Type>{
	public:
		Chain(unsigned int N, unsigned int n, unsigned int m, std::string filename);
		virtual ~Chain();

	protected:
		Vector<unsigned int> get_neighbourg(unsigned int i);
		unsigned int a_;//vector of the unit cell
};

template<typename Type>
Chain<Type>::Chain(unsigned int N, unsigned int n, unsigned int m, std::string filename):
	GenericSystem<Type>(N,n,m,2,filename),
	a_(N/m)
{
	if(a_*m != N){ std::cerr<<"Chain::Chain(N,n,m,z,filename) : The unit cell is not well defined"<<std::endl; }
	this->compute_sts();
}

template<typename Type>
Chain<Type>::~Chain(){}

template<typename Type>
Vector<unsigned int> Chain<Type>::get_neighbourg(unsigned int i){
	Vector<unsigned int> neighbourg(this->z_);
	if( i != this->n_-1){ neighbourg(0) = i+1;}
	else { neighbourg(0) = 0;}
	if( i != 0){ neighbourg(1) = i-1;}
	else { neighbourg(1) = this->n_-1;}

	return neighbourg;
}
#endif
