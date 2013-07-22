#ifndef DEF_HONEYCOMB
#define DEF_HONEYCOMB

#include "CreateSystem.hpp"

template<typename Type>
class Honeycomb: public CreateSystem<Type>{
	public:
		Honeycomb(Parseur& P);
		~Honeycomb();

	protected:
		unsigned int N_row, N_col;

		void compute_H();
};

template<typename Type>
Honeycomb<Type>::Honeycomb(Parseur& P):
	CreateSystem<Type>(P,3),
	N_row(floor(sqrt(this->N_m))),
	N_col(floor(sqrt(this->N_m)))
{
	P.set("bc",this->bc);
	if(!P.status()){
		if(this->N_m==this->N_row*this->N_col){
			this->compute_H();
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
	for(unsigned int l(0);l<this->N_row;l++){
		for(unsigned int c(0);c<this->N_col;c++){
			this->H(i,i+1) = 1;
			this->H(i,i+3) = 1;
			i++;//1
			this->H(i,i+1) = 1; 
			i++;//2
			if(c+1==N_col){
				this->H(i+1-c*4, i) = 1;
			} else {
				this->H(i,i+5) = 1;
			}
			if(l+1==N_row){
				this->H(i-1-l*N_col*4,i) = 1; 
			} else {
				this->H(i,i-1+N_col*4) = 1;
			}
			i++;//3
			if(l+1==N_row){
				this->H(i-3-l*N_col*4, i) = 1; 
			} else {
				this->H(i,i-3+N_col*4) = 1;
			}
			i++;//4
		}
	}
	this->H += this->H.transpose();
}
#endif


