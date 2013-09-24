#ifndef DEF_HONEYCOMB
#define DEF_HONEYCOMB

#include "CreateSystem.hpp"

template<typename Type>
class Honeycomb: public CreateSystem<Type>{
	public:
		Honeycomb(Parseur& P);
		~Honeycomb();

	protected:
		unsigned int Ly_, Lx_;

		void compute_H();
};

template<typename Type>
Honeycomb<Type>::Honeycomb(Parseur& P):
	CreateSystem<Type>(P,3),
	Ly_(floor(sqrt(this->m_))),
	Lx_(floor(sqrt(this->m_)))
{
	P.set("bc",this->bc_);
	if(!P.status()){
		if(this->m_==Ly_*Lx_){
			compute_H();
			this->compute_sts();
		} else {
			std::cerr<<"Honeycomb<Type> : the cluster is not a square"<<std::endl;
		}
	}
}

template<typename Type>
Honeycomb<Type>::~Honeycomb(){}

template<typename Type>
void Honeycomb<Type>::compute_H(){
	unsigned int i(0);
	for(unsigned int l(0);l<Ly_;l++){
		for(unsigned int c(0);c<Lx_;c++){
			this->H_(i,i+1) = 1;
			this->H_(i,i+3) = 1;
			i++;//1
			this->H_(i,i+1) = 1; 
			i++;//2
			if(c+1==Lx_){
				this->H_(i+1-c*4, i) = 1;
			} else {
				this->H_(i,i+5) = 1;
			}
			if(l+1==Ly_){
				this->H_(i-1-l*Lx_*4,i) = 1; 
			} else {
				this->H_(i,i-1+Lx_*4) = 1;
			}
			i++;//3
			if(l+1==Ly_){
				this->H_(i-3-l*Lx_*4, i) = 1; 
			} else {
				this->H_(i,i-3+Lx_*4) = 1;
			}
			i++;//4
		}
	}
	this->H_ += this->H_.transpose();
}
#endif


