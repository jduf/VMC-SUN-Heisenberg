#ifndef DEF_CHAIN
#define DEF_CHAIN

#include "CreateSystem.hpp"

template<typename Type>
class Chain: public CreateSystem<Type>{
	public:
		Chain(Parseur& P, std::string filename);
		virtual ~Chain();

	protected:
		Vector<unsigned int> get_neighbourg(unsigned int i);
};

template<typename Type>
Chain<Type>::Chain(Parseur& P, std::string filename):
	CreateSystem<Type>(P,2,filename)
{
	this->ref_(1) = 1;
	this->ref_(0) = 2;
	if(!P.status()){
		this->compute_sts();
		if(this->M_ % 2 == 0){ 
			this->filename_ += "-A";
			this->bc_ = -1;
		} else {
			this->filename_ += "-P";
			this->bc_ = 1;
		}
	}
}

template<typename Type>
Chain<Type>::~Chain(){ }

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
