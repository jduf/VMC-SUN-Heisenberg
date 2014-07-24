#ifndef DEF_SQUARE
#define DEF_SQUARE

#include "System2D.hpp"

template<typename Type>
class Square: public System2D<Type>{
	public:
		/*{Description*/
		/*!Constructor that organises the n sites according to the ratio Lx/Ly
		 * for a system with spuc sites per unit cell. Calls the
		 * GenericSystem<Type>(4,filename), to construct a system with 4 links
		 * per sites */
		/*}*/
		Square(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, std::string const& filename);
		virtual ~Square()=0;

	protected:
		Matrix<int> get_neighbourg(unsigned int i) const;
};

template<typename Type>
Square<Type>::Square(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(Lx,Ly,spuc,4,filename,0,0)
{
	std::cerr<<"Square::Square(N,n,m,filename) : need to set the boundary condition, and check everything"<<std::endl;
	if(this->n_==this->Ly_*this->Lx_){
		this->filename_ += "-" + tostring(this->Lx_) + "x" + tostring(this->Ly_);
		this->compute_links();
		this->status_--;
	} else {
		std::cerr<<"Square<Type> : the cluster is impossible, n must be a"<<std::endl; 
		std::cerr<<"             : multiple of "<<Lx*Ly*this->spuc_<<" ("<<Lx<<"x"<<Ly<<"x"<<this->spuc_<<")"<<std::endl; 
	}
}

template<typename Type>
Square<Type>::~Square(){}

template<typename Type>
Matrix<int> Square<Type>::get_neighbourg(unsigned int i) const {
	Matrix<int> nb(this->z_,2,1);
	/*+x neighbour*/
	if((i+1)%this->Lx_){ nb(0,0) = i+1; } 
	else { 
		nb(0,0) = (i/this->Lx_)*this->Lx_; 
		nb(0,1) = this->bc_; 
	}
	/*+y neighbour*/
	if(i<this->n_-this->Lx_){ nb(1,0) = i+this->Lx_; }
	else { 
		nb(1,0) = i-this->n_+this->Lx_; 
		nb(1,1) = this->bc_; 
	}
	/*-x neighbour*/
	if(i%this->Lx_){ nb(2,0) = i-1; }
	else {
		nb(2,0) = i+this->Lx_-1;
		nb(2,1) = this->bc_; 
	}
	/*-y neighbour*/
	if(i>=this->Lx_){ nb(3,0) = i-this->Lx_; }
	else { 
		nb(3,0) = this->n_-this->Lx_+i; 
		nb(3,1) = this->bc_; 
	}
	return nb;
}
#endif
