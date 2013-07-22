#ifndef DEF_SQUARE
#define DEF_SQUARE

#include "CreateSystem.hpp"

template<typename Type>
class Square: public CreateSystem<Type>{
	public:
		Square(Parseur& P);
		~Square();

	protected:
		unsigned int N_row, N_col;

		void compute_H();
		void save(std::string filename);
};

template<typename Type>
Square<Type>::Square(Parseur& P):
	CreateSystem<Type>(P,4),
	N_row(floor(sqrt(this->N_site))),
	N_col(floor(sqrt(this->N_site)))
{
	P.set("bc",this->bc);
	if(!P.status()){
		if(this->N_site==this->N_row*this->N_col){
			compute_H();
			this->compute_sts();
		} else {
			std::cerr<<"Square : the cluster is not a square"<<std::endl;
		}
	}
}

template<typename Type>
Square<Type>::~Square(){}

template<typename Type>
void Square<Type>::compute_H(){
	for(unsigned int i(0); i < this->N_site; i++){
		/*horizontal hopping*/
		if( (i+1) % this->N_col ){ this->H(i,i+1) = 1;}	
		else { this->H(i+1-this->N_col,i) = 1;}
		/*vertical hopping*/
		if( i+this->N_col < this->N_site ){  this->H(i,i+this->N_col) = 1; } 
		else { this->H(i-(this->N_row-1)*this->N_col,i) = 1;}
	}
	//this->H += this->H.transpose();
}
#endif

