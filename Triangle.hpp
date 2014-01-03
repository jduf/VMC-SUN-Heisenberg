#ifndef DEF_TRIANGLE
#define DEF_TRIANGLE

#include "CreateSystem.hpp"
#include "PSTricks.hpp"

template<typename Type>
class Triangle: public CreateSystem<Type>{
	public:
		Triangle(Parseur& P, std::string filename);
		~Triangle();

	protected:
		unsigned int Lx_;//!< dimension of the lattice along x-axis
		unsigned int Ly_;//!< dimension of the lattice along y-axis
		Matrix<Type> Px_;//!< translation operator along x-axis 
		Matrix<Type> Py_;//!< translation operator along y-axis 

		void compute_H();
		double occupation_number(Vector<double>& ni);
};

template<typename Type>
Triangle<Type>::Triangle(Parseur& P, std::string filename):
	CreateSystem<Type>(P,6,filename),
	Lx_(std::floor(std::sqrt(this->n_))),
	Ly_(std::floor(std::sqrt(this->n_)))
{
	this->bc_= P.get<double>("bc");
	if(!P.status()){
		if(this->n_==Ly_*Lx_){
			compute_H();
			this->compute_sts();
		} else {
			std::cerr<<"Triangle : the cluster is not a triangle"<<std::endl;
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
				if( (i+1) % Lx_ ){this->H_(i-Lx_*(Ly_-1)+1,i) = -4; }/*y jump across boundary*/
				if( i+Lx_ < this->n_ ){ this->H_(i,i+1) = -3;}/*x jump across boundary*/
			} else {
				this->H_(0,this->n_-1) = 2;
			}
		}
	}
	this->H_ += this->H_.transpose();
}

template<typename Type>
double Triangle<Type>::occupation_number(Vector<double>& ni){
	double max(0);
	for(unsigned int i(0);i<Lx_;i++){
		for(unsigned int j(0);j<Ly_;j++){
			for(unsigned int k(0);k<this->m_;k++){
				ni(i+j*Lx_) += norm_squared(this->T_(i+j*Lx_,k));
			}
			if(ni(i+j*Lx_) > max){max = ni(i+j*Lx_);}
		}
	}
	return max;
}
#endif
