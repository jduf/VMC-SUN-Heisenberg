#ifndef DEF_SQUARE
#define DEF_SQUARE

#include "GenericSystem.hpp"

template<typename Type>
class Square: public GenericSystem<Type>{
	public:
		Square(std::string const& filename);
		virtual ~Square();

	protected:
		unsigned int Lx_;//!< dimension of the lattice along x-axis
		unsigned int Ly_;//!< dimension of the lattice along y-axis

		double occupation_number(Vector<double>& ni);
		Matrix<int> get_neighbourg(unsigned int i) const;
};

template<typename Type>
Square<Type>::Square(std::string const& filename):
	GenericSystem<Type>(4,filename),
	Lx_(std::floor(std::sqrt(this->n_))),
	Ly_(std::floor(std::sqrt(this->n_)))
{
	std::cerr<<"Square::Square(N,n,m,filename) : need to set the boundary condition"<<std::endl;
	if(this->n_==Ly_*Lx_){
		this->filename_ += "-" + tostring(Lx_) + "x" + tostring(Ly_);
		this->compute_links();
	} else {
		std::cerr<<"Square::Square(N,n,m,bc,filename) : the cluster is not a square"<<std::endl;
	}
}

template<typename Type>
Square<Type>::~Square(){}

template<typename Type>
Matrix<int> Square<Type>::get_neighbourg(unsigned int i) const {
	Matrix<int> nb(this->z_,2,1);
	/*+x neighbour*/
	if((i+1)%Lx_){ nb(0,0) = i+1; } 
	else { 
		nb(0,0) = (i/Lx_)*Lx_; 
		nb(0,1) = this->bc_; 
	}
	/*+y neighbour*/
	if(i<this->n_-Lx_){ nb(1,0) = i+Lx_; }
	else { 
		nb(1,0) = i-this->n_+Lx_; 
		nb(1,1) = this->bc_; 
	}
	/*-x neighbour*/
	if(i%Lx_){ nb(2,0) = i-1; }
	else {
		nb(2,0) = i+Lx_-1;
		nb(2,1) = this->bc_; 
	}
	/*-y neighbour*/
	if(i>=Lx_){ nb(3,0) = i-Lx_; }
	else { 
		nb(3,0) = this->n_-Lx_+i; 
		nb(3,1) = this->bc_; 
	}
	return nb;
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
