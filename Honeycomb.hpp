#ifndef DEF_HONEYCOMB
#define DEF_HONEYCOMB

#include "CreateSystem.hpp"

template<typename Type>
class Honeycomb: public CreateSystem<Type>{
	public:
		Honeycomb(Parseur& P, std::string filename);
		~Honeycomb();

	protected:
		unsigned int Lx_;//!< dimension of the lattice along x-axis
		unsigned int Ly_;//!< dimension of the lattice along y-axis
		Matrix<Type> Px_;//!< translation operator along x-axis 
		Matrix<Type> Py_;//!< translation operator along y-axis 

		void compute_H();
};

template<typename Type>
Honeycomb<Type>::Honeycomb(Parseur& P, std::string filename):
	CreateSystem<Type>(P,3,filename),
	Lx_(floor(sqrt(this->m_))),
	Ly_(floor(sqrt(this->m_)))
{
	this->bc_= P.get<double>("bc");
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
			//0
			this->H_(i,i+1) = 1;
			if(l+1<Ly_){
				this->H_(i,i+1+Lx_*4) = 1;
			} else {
				this->H_(i+1-l*Lx_*4,i) = 1;
			}
			if(c==0){
				this->H_(i,i+Lx_*4-1) = 1;
			} else {
				this->H_(i-1,i) = 1;
			}
			i+=2;//2
			this->H_(i,i-1) = 1;
			this->H_(i,i+1) = 1; 
			if(l==0){
				this->H_(i,i+1+(Ly_-1)*Lx_*4) = 1;
			} else {
				this->H_(i,i+1-Lx_*4) = 1;
			}
			i+=2;//4
		}
	}
	this->H_ += this->H_.transpose();
}
#endif


