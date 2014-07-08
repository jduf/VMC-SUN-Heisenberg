#ifndef DEF_TRIANGLE
#define DEF_TRIANGLE

#include "GenericSystem.hpp"

template<typename Type>
class Triangle: public GenericSystem<Type>{
	public:
		/*{Description*/
		/*!Constructor that organises the n sites according to the ratio Lx/Ly
		 * for a system with spuc sites per unit cell. Calls the
		 * GenericSystem<Type>(6,filename), to construct a system with 6 links
		 * per sites */
		/*}*/
		Triangle(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, std::string const& filename);
		virtual ~Triangle()=0;

	protected:
		unsigned int Lx_;//!< dimension of the lattice along x-axis
		unsigned int Ly_;//!< dimension of the lattice along y-axis
		unsigned int spuc_;//!< site per unit cell

		double occupation_number(Vector<double>& ni);
		Matrix<int> get_neighbourg(unsigned int i) const;
};

template<typename Type>
Triangle<Type>::Triangle(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, std::string const& filename):
	GenericSystem<Type>(6,filename),
	Lx_(sqrt(Lx*this->n_/(Ly*spuc))),
	Ly_(sqrt(Ly*this->n_/(Lx*spuc))),
	spuc_(spuc)
{
	std::cerr<<"Triangle::Triangle(N,n,m,filename) : need to set the boundary condition"<<std::endl;
	std::cerr<<"Triangle::Triangle(N,n,m,filename) : and need to check that they are correct"<<std::endl;
	if(this->n_==Ly_*Lx_){
		this->filename_ += "-" + tostring(Lx_) + "x" + tostring(Ly_);
		this->compute_links();
		this->status_--;
	} else {
		std::cerr<<"Triangle<Type> : the cluster is impossible, n must be a"<<std::endl; 
		std::cerr<<"               : multiple of "<<Lx*Ly*spuc_<<" ("<<Lx<<"x"<<Ly<<"x"<<spuc_<<")"<<std::endl; 
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
