#ifndef DEF_KAGOME
#define DEF_KAGOME

#include "GenericSystem.hpp"

template<typename Type>
class Kagome: public GenericSystem<Type>{
	public:
		Kagome(
				unsigned int const& N,
				unsigned int const& n, 
				unsigned int const& m, 
				int const& bc, 
				Vector<unsigned int> const& ref,  
				unsigned int const& Lx,
				unsigned int const& Ly,
				unsigned int const& spuc, 
				std::string const& filename);
		virtual ~Kagome();

	protected:
		unsigned int Lx_;//!< dimension of the lattice along x-axis
		unsigned int Ly_;//!< dimension of the lattice along y-axis
		unsigned int spuc_;

		virtual Matrix<int> get_neighbourg(unsigned int i) const;
};

template<typename Type>
Kagome<Type>::Kagome(
		unsigned int const& N, 
		unsigned int const& n, 
		unsigned int const& m, 
		int const& bc, 
		Vector<unsigned int> const& ref,  
		unsigned int const& Lx,
		unsigned int const& Ly,
		unsigned int const& spuc,
		std::string const& filename):
	GenericSystem<Type>(N,n,m,bc,ref,4,filename),
	Lx_(Lx),
	Ly_(Ly),
	spuc_(spuc)
{
	unsigned int nc(sqrt(this->n_/spuc_));
	Lx_ *= nc;
	Ly_ *= nc;
	if(this->n_==Ly_*Lx_*spuc_){
		this->filename_ += "-" + tostring(Lx_) +"x"+ tostring(Ly_);
		this->compute_links();
	} else {
		std::cerr<<"Kagome<Type> : the cluster not possible"<<std::endl;
	}
}

template<typename Type>
Kagome<Type>::~Kagome(){}

template<typename Type>
Matrix<int> Kagome<Type>::get_neighbourg(unsigned int i) const {
	Matrix<int> nb(this->z_,2,1);
	switch(spuc_){
		case 3:
			{
				switch(i%3){
					case 0:
						{
							/*+x neighbour*/
							nb(0,0) = i+1;
							/*+x+y neighbour*/
							nb(1,0) = i+2;
						}break;
					case 1:
						{
							/*+x neighbour*/
							if((i+2)%(3*Lx_)){ nb(0,0) = i+2; }
							else {
								nb(0,0) = i+2-3*Lx_;
								nb(0,1) = this->bc_;
							}
							/*-x+y neighbour*/
							nb(1,0) = i+1;
						}break;
					case 2:
						{
							/*+x+y neighbour*/
							if(i<this->n_-3*Lx_){ nb(0,0) = i+3*Lx_-2; }
							else { 
								nb(0,0) = i-2-3*Lx_*(Ly_-1);
								nb(0,1) = this->bc_;
							}
							/*-x+y neighbour*/
							if((i-2)%(3*Lx_)==0){
								if(i<this->n_-3*Lx_){
									nb(1,0) = i-4+6*Lx_;
									nb(1,1) = this->bc_;
								} else {
									nb(1,0) = 3*Lx_-2;
									nb(1,1) = this->bc_*this->bc_;
								}
							} else {
								if(i<this->n_-3*Lx_){ nb(1,0) = i-4+3*Lx_; }
								else {
									nb(1,0) = i-4-3*Lx_*(Ly_-1);
									nb(1,1) = this->bc_;
								}
							}
						}break;
				}
			}break;
		case 9:
			{
				switch(i%9){
					case 0:
						{ 
							nb(0,0) = i+1; 
							if(i%(Lx_*9)==0){ nb(1,0) = i-1+Lx_*9; }
							else { 
								nb(1,0) = i-1; 
								nb(1,1) = this->bc_; 
							}
							if(i>=(Lx_*9)){ nb(2,0) = i+6-Lx_*9; }
							else { 
								nb(2,0) = i+6+(Ly_-1)*Lx_*9;
								nb(2,1) = this->bc_; 
							}
							nb(3,0) = i+5;
						}break;
					case 1:
						{
							nb(0,0) = i+1;
							if((i-1)%(Lx_*9)!=0){ nb(1,0) = i-3; }
							else {
								nb(1,0) = i-3+Lx_*9; 
								nb(1,1) = this->bc_; 
							}
							if((i-1)%(Lx_*9)==0){ nb(2,0) = i-2+Lx_*9; }
							else { 
								nb(2,0) = i-2; 
								nb(2,1) = this->bc_; 
							}
							nb(3,0) = i-1;
						}break;
					case 2:
						{
							nb(0,0) = i+1;
							nb(1,0) = i+4;
							if((i-2)%(Lx_*9)!=0){ nb(2,0) = i-4; }
							else {
								nb(2,0) = i-4+Lx_*9; 
								nb(2,1) = this->bc_; 
							}
							nb(3,0) = i-1;
						}break;
					case 3:
						{
							nb(0,0) = i+1;
							nb(1,0) = i+5;
							nb(2,0) = i+3;
							nb(3,0) = i-1;
						}break;
					case 4:
						{
							nb(0,0) = i+1;
							if(i>=(Lx_*9)){ nb(1,0) = i+3-Lx_*9; }
							else { 
								nb(1,0) = i+3+(Ly_-1)*Lx_*9;
								nb(1,1) = this->bc_; 
							}
							nb(2,0) = i+4;
							nb(3,0) = i-1;
						}break;
					case 5:
						{
							nb(0,0) = i-5;
							if(i>=(Lx_*9)){ nb(2,0) = i+1-Lx_*9; }
							else { 
								nb(2,0) = i+1+(Ly_-1)*Lx_*9;
								nb(2,1) = this->bc_; 
							}
							if(i>=(Lx_*9)){ nb(3,0) = i+2-Lx_*9; }
							else { 
								nb(3,0) = i+2+(Ly_-1)*Lx_*9;
								nb(3,1) = this->bc_; 
							}
							nb(3,0) = i-1;
						}break;
					case 6:
						{
							if(i<(Ly_-1)*Lx_*9){ nb(0,0) = i-1+Lx_*9; }
							else {
								nb(0,0) = i-1-(Ly_-1)*Lx_*9;
								nb(0,1) = this->bc_; 
							}
							if(i<(Ly_-1)*Lx_*9){ nb(1,0) = i-6+Lx_*9; }
							else {
								nb(1,0) = i-6-(Ly_-1)*Lx_*9;
								nb(1,1) = this->bc_; 
							}
							nb(2,0) = i-4;
							nb(3,0) = i-3;
						}break;
					case 7:
						{
							if((i+2)%(Lx_*9)!=0){ nb(0,0) = i+3; }
							else {
								nb(0,0) = i+3-Lx_*9;
								nb(0,1) = this->bc_; 
							}
							if((i+2)%(Lx_*9)!=0){ nb(1,0) = i+4; }
							else {
								nb(1,0) = i+4-Lx_*9;
								nb(1,1) = this->bc_; 
							}
							if(i<(Ly_-1)*Lx_*9){ nb(2,0) = i-3+Lx_*9; }
							else {
								nb(2,0) = i-3-(Ly_-1)*Lx_*9;
								nb(2,1) = this->bc_; 
							}
							if(i<(Ly_-1)*Lx_*9){ nb(3,0) = i-2+Lx_*9; }
							else {
								nb(3,0) = i-2-(Ly_-1)*Lx_*9;
								nb(3,1) = this->bc_; 
							}
						}break;
					case 8:
						{
							if((i+1)%(Lx_*9)!=0){ nb(0,0) = i+1; }
							else {
								nb(0,0) = i+1-Lx_*9;
								nb(0,1) = this->bc_; 
							}
							if((i+1)%(Lx_*9)!=0){ nb(1,0) = i+2; }
							else {
								nb(1,0) = i+2-Lx_*9;
								nb(1,1) = this->bc_; 
							}
							nb(2,0) = i-5;
							nb(3,0) = i-4;
						}break;
				}
			}break;
	}
	return nb;
}
#endif
