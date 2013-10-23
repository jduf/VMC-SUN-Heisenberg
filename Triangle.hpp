#ifndef DEF_TRIANGLE
#define DEF_TRIANGLE

#include "CreateSystem.hpp"

template<typename Type>
class Triangle: public CreateSystem<Type>{
	public:
		Triangle(Parseur& P);
		~Triangle();

	protected:
		unsigned int Lx_;//!< dimension of the lattice along x-axis
		unsigned int Ly_;//!< dimension of the lattice along y-axis
		Matrix<Type> Px_;//!< translation operator along x-axis 
		Matrix<Type> Py_;//!< translation operator along y-axis 

		void compute_H();
		void save(std::string filename);
};

template<typename Type>
Triangle<Type>::Triangle(Parseur& P):
	CreateSystem<Type>(P,6),
	Lx_(std::floor(std::sqrt(this->n_))),
	Ly_(std::floor(std::sqrt(this->n_)))
{
	this->bc_= P.get<double>("bc");
	if(!P.status()){
		if(this->n_==Ly_*Lx_){
			compute_H();
			this->compute_sts();
		} else {
			std::cerr<<"Triangle : the cluster is not a Triangle"<<std::endl;
		}
	}
}

template<typename Type>
Triangle<Type>::~Triangle(){}

template<typename Type>
void Triangle<Type>::compute_H(){
	for(unsigned int i(0); i < this->n_; i++){
		/*horizontal hopping*/
		if( (i+1) % Lx_ ){ this->H_(i,i+1) = 1;}	
		else { this->H_(i+1-Lx_,i) = -1;}
		/*vertical hopping*/
		if( i+Lx_ < this->n_ ){  this->H_(i,i+Lx_) = 1; } 
		else { this->H_(i-(Ly_-1)*Lx_,i) = -2;}
		/*diagonal hopping*/
		if( (i+1) % Lx_ && i+Lx_ < this->n_ ){  this->H_(i,i+Lx_+1) = 1; } 
		else {
			if(i+1 < this->n_ ){
				if( !((i+1) % Lx_) ){ this->H_(i,i+1) = -1;}/*x jump across boundary*/
				if( i+Lx_ >= this->n_ ){  this->H_(i-Lx_*(Ly_-1)+1,i) = 1; }/*y jump across boundary*/
			} else {
				this->H_(0,this->n_-1) = -3;
			}
		}
	}
	this->H_ += this->H_.transpose();
}

#endif

