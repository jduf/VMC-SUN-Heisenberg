#ifndef DEF_TRIANGLE
#define DEF_TRIANGLE

#include "System2D.hpp"

template<typename Type>
class Triangle: public System2D<Type>{
	public:
		/*{Description*/
		/*!Constructor that organises the n sites according to the ratio Lx/Ly
		 * for a system with spuc sites per unit cell. Calls the
		 * GenericSystem<Type>(6,filename), to construct a system with 6 links
		 * per sites */
		/*}*/
		Triangle(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Triangle()=0;

	protected:
		/*!Returns the neighbours of site i*/
		Matrix<int> get_neighbourg(unsigned int i) const;
};

template<typename Type>
Triangle<Type>::Triangle(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(Lx,Ly,spuc,6,filename,0,0)
{
	std::cerr<<"Triangle::Triangle(N,n,m,filename) : need to set the boundary condition"<<std::endl;
	std::cerr<<"Triangle::Triangle(N,n,m,filename) : and need to check that they are correct"<<std::endl;
	if(this->status_==2){ this->compute_links(); }
}

template<typename Type>
Triangle<Type>::~Triangle(){}

template<typename Type>
Matrix<int> Triangle<Type>::get_neighbourg(unsigned int i) const {
	Matrix<int> nb(this->z_,2,1);
	/*0 neighbour*/
	if((i+1)%this->Lx_){ nb(0,0) = i+1; } 
	else {
		nb(0,0) = (i/this->Lx_)*this->Lx_;
		nb(0,1) = this->bc_; 
	}
	/*pi/3 neighbour*/
	if((i+1)%this->Lx_ && i+this->Lx_<this->n_){ nb(1,0) = i+this->Lx_+1; } 
	else {
		if(i+1<this->n_){
			if((i+1)%this->Lx_){
				nb(1,0) = i-this->Lx_*(this->Ly_-1)+1; 
			}
			if(i+this->Lx_<this->n_ ){ 
				nb(1,0) = i+1; 
				nb(1,1) = this->bc_; 
			}
		} else {
			nb(1,0) = 0 ;
			nb(1,1) = this->bc_; 
		}
	}
	/*2pi/3 neighbour*/
	if(i+this->Lx_<this->n_){ nb(2,0) = i+this->Lx_; }
	else {
		nb(2,0) = i-this->n_+this->Lx_;
		nb(2,1) = this->bc_; 
	}
	/*pi neighbour*/
	if(i%this->Lx_){ nb(3,0) = i-1; }
	else {
		nb(3,0) = i+this->Lx_-1;
		nb(3,1) = this->bc_; 
	}
	/*4pi/3 neighbour*/
	if(i%this->Lx_ && i>=this->Lx_){ nb(4,0) = i-this->Lx_-1; } 
	else {
		if(i!=0){
			if(i<this->Lx_){
				nb(4,0) = this->n_-this->Lx_+i-1;
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
	if(i>=this->Lx_){ nb(5,0) = i-this->Lx_; }
	else {
		nb(5,0) = this->n_-this->Lx_+i; 
		nb(5,1) = this->bc_; 
	}

	return nb;
}
#endif
