#ifndef DEF_HONEYCOMB
#define DEF_HONEYCOMB

#include "System2D.hpp"

template<typename Type>
class Honeycomb: public System2D<Type>{
	public:
		/*{Description*/
		/*!Constructor that organises the n sites according to the ratio Lx/Ly
		 * for a system with spuc sites per unit cell. Calls the
		 * GenericSystem<Type>(3,filename), to construct a system with 3 links
		 * per sites */
		/*}*/
		Honeycomb(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, std::string const& filename);
		/*!Pure virtual destructor (abstract class)*/
		virtual ~Honeycomb()=0;

	protected:
		/*!Returns the neighbours of site i*/
		Matrix<int> get_neighbourg(unsigned int i) const;
};

template<typename Type>
Honeycomb<Type>::Honeycomb(unsigned int const& Lx, unsigned int const& Ly, unsigned int const& spuc, std::string const& filename):
	System2D<Type>(Lx,Ly,spuc,3,filename)
{
	if(this->status_==2){ this->compute_links(); }
}

template<typename Type>
Honeycomb<Type>::~Honeycomb(){}

template<typename Type>
Matrix<int> Honeycomb<Type>::get_neighbourg(unsigned int i) const {
	Matrix<int> nb(this->z_,2,1);
	switch(this->spuc_){
		case 6:
			{
				switch(i%6){
					case 0:
						{
							/*+x+y neighbour*/
							nb(0,0) = i+1;
							/*-x neighbour*/
							if(i%(this->spuc_*this->Lx_)){ nb(1,0) = i-3; }
							else {
								nb(1,0) = i-3+this->spuc_*this->Lx_;
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
							if(i<this->n_-this->spuc_*this->Lx_){ 
								if((i-1)%(this->spuc_*this->Lx_)){ nb(1,0) = i-3+this->spuc_*this->Lx_ ; }
								else {
									nb(1,0) = i-3+2*this->spuc_*this->Lx_;
									nb(1,1) = this->bc_;
								}
							} else {
								if(i-1!=this->n_-this->spuc_*this->Lx_){
									nb(1,0) = i-3-this->spuc_*this->Lx_*(this->Ly_-1);
									nb(1,1) = this->bc_;
								} else {
									nb(1,0) = -2+this->spuc_*this->Lx_;
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
							if(i<this->n_-this->spuc_*this->Lx_){ nb(1,0) = i+3+this->spuc_*this->Lx_; }
							else {
								nb(1,0) = i+3-this->spuc_*this->Lx_*(this->Ly_-1);
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
							if((i+3)%(this->spuc_*this->Lx_)){ nb(1,0) = i+3; }
							else {
								nb(1,0) = i+3-this->spuc_*this->Lx_;
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
							if(i>this->spuc_*this->Lx_){ 
								if((i+2)%(this->spuc_*this->Lx_)){ nb(1,0) = i+3-this->spuc_*this->Lx_; }
								else {
									nb(1,0) = i+3-2*this->spuc_*this->Lx_;
									nb(1,1) = this->bc_;
								}
							} else {
								if(i+2!=this->spuc_*this->Lx_){
									nb(1,0) = i+3+this->spuc_*this->Lx_*(this->Ly_-1);
									nb(1,1) = this->bc_;
								} else {
									nb(1,0) = this->n_+1-this->spuc_*this->Lx_;
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
							if(i>this->spuc_*this->Lx_){ nb(1,0) = i-3-this->spuc_*this->Lx_; }
							else {
								nb(1,0) = i-3+this->spuc_*this->Lx_*(this->Ly_-1);
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
