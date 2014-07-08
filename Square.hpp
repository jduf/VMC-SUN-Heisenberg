#ifndef DEF_SQUARE
#define DEF_SQUARE

#include "GenericSystem.hpp"

template<typename Type>
class Square: public GenericSystem<Type>{
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
		unsigned int Lx_;//!< dimension of the lattice along x-axis
		unsigned int Ly_;//!< dimension of the lattice along y-axis
		unsigned int spuc_;//!< site per unit cell

		Matrix<int> get_neighbourg(unsigned int i) const;
		double occupation_number(Vector<double>& ni);
};

template<typename Type>
Square<Type>::Square(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, std::string const& filename):
	GenericSystem<Type>(4,filename),
	Lx_(sqrt(Lx*this->n_/(Ly*spuc))),
	Ly_(sqrt(Ly*this->n_/(Lx*spuc))),
	spuc_(spuc)
{
	std::cerr<<"Square::Square(N,n,m,filename) : need to set the boundary condition, and check everything"<<std::endl;
	if(this->n_==Ly_*Lx_){
		this->filename_ += "-" + tostring(Lx_) + "x" + tostring(Ly_);
		this->compute_links();
		this->status_--;
	} else {
		std::cerr<<"Square<Type> : the cluster is impossible, n must be a"<<std::endl; 
		std::cerr<<"             : multiple of "<<Lx*Ly*spuc_<<" ("<<Lx<<"x"<<Ly<<"x"<<spuc_<<")"<<std::endl; 
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
