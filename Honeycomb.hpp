#ifndef DEF_HONEYCOMB
#define DEF_HONEYCOMB

#include "GenericSystem.hpp"

template<typename Type>
class Honeycomb: public GenericSystem<Type>{
	public:
		/*{Description*/
		/*!Constructor that organises the n sites according to the ratio Lx/Ly
		 * for a system with spuc sites per unit cell. Calls the
		 * GenericSystem<Type>(3,filename), to construct a system with 3 links
		 * per sites */
		/*}*/
		Honeycomb(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, std::string const& filename);
		virtual ~Honeycomb()=0;

	protected:
		unsigned int Lx_;	//!< dimension of the lattice along x-axis
		unsigned int Ly_;	//!< dimension of the lattice along y-axis
		unsigned int spuc_;	//!< site per unit cell

		Matrix<int> get_neighbourg(unsigned int i) const;
};

template<typename Type>
Honeycomb<Type>::Honeycomb(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, std::string const& filename):
	GenericSystem<Type>(3,filename),
	Lx_(sqrt(Lx*this->n_/(Ly*spuc))),
	Ly_(sqrt(Ly*this->n_/(Lx*spuc))),
	spuc_(spuc)
{
	if(this->n_==Ly_*Lx_*spuc_){
		this->filename_ += "-" + tostring(Lx_) +"x"+ tostring(Ly_);
		this->compute_links();
		this->status_--;
	} else {
		std::cerr<<"Honeycomb<Type> : the cluster is impossible, n must be a"<<std::endl; 
		std::cerr<<"                : multiple of "<<Lx*Ly*spuc_<<" ("<<Lx<<"x"<<Ly<<"x"<<spuc_<<")"<<std::endl; 
	}
}

template<typename Type>
Honeycomb<Type>::~Honeycomb(){}

template<typename Type>
Matrix<int> Honeycomb<Type>::get_neighbourg(unsigned int i) const {
	Matrix<int> nb(this->z_,2,1);
	switch(spuc_){
		case 6:
			{
				switch(i%6){
					case 0:
						{
							/*+x+y neighbour*/
							nb(0,0) = i+1;
							/*-x neighbour*/
							if(i%(spuc_*Lx_)){ nb(1,0) = i-3; }
							else {
								nb(1,0) = i-3+spuc_*Lx_;
								nb(1,1) = this->bc_;
							}
							/*+x-y neighbour*/
							nb(2,0) = i+5;
						}break;
					case 1:
						{
							/*+x neighbour*/
							nb(0,0) = i+1;
							/*-x+y neighbour*/
							if(i<this->n_-spuc_*Lx_){ 
								if((i-1)%(spuc_*Lx_)){ nb(1,0) = i-3+spuc_*Lx_ ; }
								else {
									nb(1,0) = i-3+2*spuc_*Lx_;
									nb(1,1) = this->bc_;
								}
							} else {
								if(i-1!=this->n_-spuc_*Lx_){
									nb(1,0) = i-3-spuc_*Lx_*(Ly_-1);
									nb(1,1) = this->bc_;
								} else {
									nb(1,0) = -2+spuc_*Lx_;
									nb(1,1) = this->bc_*this->bc_;
								}
							}
							/*-x-y neighbour*/
							nb(2,0) = i-1;
						}break;
					case 2:
						{
							/*+x-y neighbour*/
							nb(0,0) = i+1;
							/*+x+y neighbour*/
							if(i<this->n_-spuc_*Lx_){ nb(1,0) = i+3+spuc_*Lx_; }
							else {
								nb(1,0) = i+3-spuc_*Lx_*(Ly_-1);
								nb(1,1) = this->bc_;
							}
							/*-x neighbour*/
							nb(2,0) = i-1;
						}break;
					case 3:
						{
							/*-x-y neighbour*/
							nb(0,0) = i+1;
							/*+x neighbour*/
							if((i+3)%(spuc_*Lx_)){ nb(1,0) = i+3; }
							else {
								nb(1,0) = i+3-spuc_*Lx_;
								nb(1,1) = this->bc_;
							}
							/*-x+y neighbour*/
							nb(2,0) = i-1;
						}break;
					case 4:
						{
							/*-x neighbour*/
							nb(0,0) = i+1;
							/*+x-y neighbour*/
							if(i>spuc_*Lx_){ 
								if((i+2)%(spuc_*Lx_)){ nb(1,0) = i+3-spuc_*Lx_; }
								else {
									nb(1,0) = i+3-2*spuc_*Lx_;
									nb(1,1) = this->bc_;
								}
							} else {
								if(i+2!=spuc_*Lx_){
									nb(1,0) = i+3+spuc_*Lx_*(Ly_-1);
									nb(1,1) = this->bc_;
								} else {
									nb(1,0) = this->n_+1-spuc_*Lx_;
									nb(1,1) = this->bc_*this->bc_;
								}
							}
							/*+x+y neighbour*/
							nb(2,0) = i-1;
						}break;
					case 5:
						{
							/*-x+y neighbour*/
							nb(0,0) = i-5;
							/*-x-y neighbour*/
							if(i>spuc_*Lx_){ nb(1,0) = i-3-spuc_*Lx_; }
							else {
								nb(1,0) = i-3+spuc_*Lx_*(Ly_-1);
								nb(1,1) = this->bc_;
							}
							/*+x neighbour*/
							nb(2,0) = i-1;
						}break;
				}
			}break;
	}
	return nb;
}
#endif
