#ifndef DEF_SQUARE
#define DEF_SQUARE

#include "GenericSystem.hpp"

template<typename Type>
class Square: public GenericSystem<Type>{
	public:
		Square(unsigned int N, unsigned int n, unsigned int m, std::string filename);
		virtual ~Square();

	protected:
		unsigned int Lx_;//!< dimension of the lattice along x-axis
		unsigned int Ly_;//!< dimension of the lattice along y-axis
		Matrix<unsigned int> BC_;

		double occupation_number(Vector<double>& ni);
		Vector<unsigned int> get_neighbourg(unsigned int i);
};

template<typename Type>
Square<Type>::Square(unsigned int N, unsigned int n, unsigned int m, std::string filename):
	GenericSystem<Type>(N,n,m,4,filename),
	Lx_(std::floor(std::sqrt(this->n_))),
	Ly_(std::floor(std::sqrt(this->n_)))
{
	this->bc_ = -1;
	std::cerr<<"Square::Square(N,n,m,filename) : need to set the boundary condition"<<std::endl;
	if(this->n_==Ly_*Lx_){
		this->compute_sts();
		this->filename_ += "-" + tostring(Lx_) + "x" + tostring(Ly_);
		if(this->bc_ != 0){
			if(this->bc_ == 1){ this->filename_ += "-P";} 
			else{ this->filename_ += "-A";}
			unsigned int k(0);
			BC_.set(Lx_+Ly_,2);
			for(unsigned int i(0); i < this->n_; i++){
				if(!((i+1) % Lx_ )){
					BC_(k,0) = i;
					BC_(k,1) = i+1-Lx_;
					k++;
				}	
				if( i+Lx_>=this->n_){ 
					BC_(k,0) = i;
					BC_(k,1) = i-(Ly_-1)*Lx_;
					k++;
				}
			}
		}
	} else {
		std::cerr<<"Square : the cluster is not a square"<<std::endl;
	}
}

template<typename Type>
Square<Type>::~Square(){}

template<typename Type>
Vector<unsigned int> Square<Type>::get_neighbourg(unsigned int i){
	Vector<unsigned int> neighbourg(this->z_);
	/*+x neighbour*/
	if((i+1)%Lx_){ neighbourg(0) = i+1; } 
	else { neighbourg(0) = (i/Lx_)*Lx_; }
	/*+y neighbour*/
	if(i<this->n_-Lx_){ neighbourg(1) = i+Lx_; }
	else { neighbourg(1) = i-this->n_+Lx_; }
	/*-x neighbour*/
	if(i%Lx_){ neighbourg(2) = i-1; }
	else { neighbourg(2) = i+Lx_-1; }
	/*-y neighbour*/
	if(i>=Lx_){ neighbourg(3) = i-Lx_; }
	else { neighbourg(3) = this->n_-Lx_+i; }

	return neighbourg;
}

template<typename Type>
double Square<Type>::occupation_number(Vector<double>& ni){
	double max(0);
	for(unsigned int i(0);i<Lx_;i++){
		for(unsigned int j(0);j<Ly_;j++){
			for(unsigned int k(0);k<this->M_;k++){
				ni(i+j*Lx_) += norm_squared(this->T_(i+j*Lx_,k));
			}
			if(ni(i+j*Lx_) > max){max = ni(i+j*Lx_);}
		}
	}
	return max;
}
#endif
