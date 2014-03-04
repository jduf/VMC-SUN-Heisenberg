#ifndef DEF_TRIANGLE
#define DEF_TRIANGLE

#include "GenericSystem.hpp"

template<typename Type>
class Triangle: public GenericSystem<Type>{
	public:
		Triangle(unsigned int N, unsigned int n, unsigned int m, std::string filename);
		virtual ~Triangle();

	protected:
		unsigned int Lx_;//!< dimension of the lattice along x-axis
		unsigned int Ly_;//!< dimension of the lattice along y-axis

		double occupation_number(Vector<double>& ni);
		Matrix<int> get_neighbourg(unsigned int i) const;
};

template<typename Type>
Triangle<Type>::Triangle(unsigned int N, unsigned int n, unsigned int m, std::string filename):
	GenericSystem<Type>(N,n,m,6,filename),
	Lx_(std::floor(std::sqrt(this->n_))),
	Ly_(std::floor(std::sqrt(this->n_)))
{
	this->bc_ = -1;
	std::cerr<<"Triangle::Triangle(N,n,m,filename) : need to set the boundary condition"<<std::endl;
	std::cerr<<"Triangle::Triangle(N,n,m,filename) : and need to check that they are correct"<<std::endl;
	if(this->n_==Ly_*Lx_){
		this->compute_sts();
		this->filename_ += "-" + tostring(Lx_) + "x" + tostring(Ly_);
		if(this->bc_ != 0){
			if(this->bc_ == 1){ this->filename_ += "-P";} 
			else{ this->filename_ += "-A";}
		}
	} else {
		std::cerr<<"Triangle : the cluster is not a square"<<std::endl;
	}
}

template<typename Type>
Triangle<Type>::~Triangle(){}

template<typename Type>
Matrix<int> Triangle<Type>::get_neighbourg(unsigned int i) const {
	Matrix<int> nb(this->z_,2,1);
	/*0 neighbour*/
	if((i+1)%Lx_){ nb(0,0) = i+1; } 
	else {
		nb(0,0) = (i/Lx_)*Lx_;
		nb(0,1) = this->bc_; 
	}
	/*pi/3 neighbour*/
	if((i+1)%Lx_ && i+Lx_<this->n_){ nb(1,0) = i+Lx_+1; } 
	else {
		if(i+1<this->n_){
			if((i+1)%Lx_){
				nb(1,0) = i-Lx_*(Ly_-1)+1; 
			}
			if(i+Lx_<this->n_ ){ 
				nb(1,0) = i+1; 
				nb(1,1) = this->bc_; 
			}
		} else {
			nb(1,0) = 0 ;
			nb(1,1) = this->bc_; 
		}
	}
	/*2pi/3 neighbour*/
	if(i+Lx_<this->n_){ nb(2,0) = i+Lx_; }
	else {
		nb(2,0) = i-this->n_+Lx_;
		nb(2,1) = this->bc_; 
	}
	/*pi neighbour*/
	if(i%Lx_){ nb(3,0) = i-1; }
	else {
		nb(3,0) = i+Lx_-1;
		nb(3,1) = this->bc_; 
	}
	/*4pi/3 neighbour*/
	if(i%Lx_ && i>=Lx_){ nb(4,0) = i-Lx_-1; } 
	else {
		if(i!=0){
			if(i<Lx_){
				nb(4,0) = this->n_-Lx_+i-1;
				nb(4,1) = this->bc_; 
			} else{
				nb(4,0) = i-1;
				nb(4,1) = this->bc_; 
			}
		} else { 
			nb(4,0) = this->n_-1 ;
			nb(4,1) = this->bc_; 
		}
	}
	/*5pi/3 neighbour*/
	if(i>=Lx_){ nb(5,0) = i-Lx_; }
	else {
		nb(5,0) = this->n_-Lx_+i; 
		nb(5,1) = this->bc_; 
	}

	return nb;
}

template<typename Type>
double Triangle<Type>::occupation_number(Vector<double>& ni){
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
