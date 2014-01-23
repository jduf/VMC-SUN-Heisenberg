#ifndef DEF_TRIANGLE
#define DEF_TRIANGLE

#include "CreateSystem.hpp"
#include "PSTricks.hpp"

template<typename Type>
class Triangle: public CreateSystem<Type>{
	public:
		Triangle(Parseur& P, std::string filename);
		virtual ~Triangle();

	protected:
		unsigned int Lx_;//!< dimension of the lattice along x-axis
		unsigned int Ly_;//!< dimension of the lattice along y-axis
		Matrix<Type> Px_;//!< translation operator along x-axis 
		Matrix<Type> Py_;//!< translation operator along y-axis 
		Matrix<unsigned int> BC_;

		double occupation_number(Vector<double>& ni);
		Vector<unsigned int> get_neighbourg(unsigned int i);
};

template<typename Type>
Triangle<Type>::Triangle(Parseur& P, std::string filename):
	CreateSystem<Type>(P,6,filename),
	Lx_(std::floor(std::sqrt(this->n_))),
	Ly_(std::floor(std::sqrt(this->n_))),
	BC_(Lx_+Ly_,2)
{
	this->ref_(0) = 3;
	this->bc_= P.get<double>("bc");
	if(!P.status()){
		if(this->n_==Ly_*Lx_){
			this->compute_sts();
			unsigned int k(0);
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
			this->filename_ += "-N" + tostring(this->N_);
			this->filename_ += "-S" + tostring(this->n_);
			this->filename_ += "-" + tostring(Lx_) + "x" + tostring(Ly_);
		} else {
			std::cerr<<"Triangle : the cluster is not a square"<<std::endl;
		}
	}
}

template<typename Type>
Triangle<Type>::~Triangle(){}

template<typename Type>
Vector<unsigned int> Triangle<Type>::get_neighbourg(unsigned int i){
	Vector<unsigned int> neighbourg(this->z_);
	/*0 neighbour*/
	if((i+1)%Lx_){ neighbourg(0) = i+1; } 
	else { neighbourg(0) = (i/Lx_)*Lx_; }
	/*pi/3 neighbour*/
	if((i+1)%Lx_ && i+Lx_<this->n_){ neighbourg(1) = i+Lx_+1; } 
	else {
		if(i+1<this->n_){
			if((i+1)%Lx_){ neighbourg(1) = i-Lx_*(Ly_-1)+1; }
			if(i+Lx_<this->n_ ){ neighbourg(1) = i+1; }
		} else { neighbourg(1) = 0 ;}
	}
	/*2pi/3 neighbour*/
	if(i+Lx_<this->n_){ neighbourg(2) = i+Lx_; }
	else { neighbourg(2) = i-this->n_+Lx_; }
	/*pi neighbour*/
	if(i%Lx_){ neighbourg(3) = i-1; }
	else { neighbourg(3) = i+Lx_-1; }
	/*4pi/3 neighbour*/
	if(i%Lx_ && i>=Lx_){ neighbourg(4) = i-Lx_-1; } 
	else {
		if(i!=0){
			if(i<Lx_){ neighbourg(4) = this->n_-Lx_+i-1;}
			else{ neighbourg(4) = i-1; }
		} else { neighbourg(4) = this->n_-1 ;}
	}
	/*5pi/3 neighbour*/
	if(i>=Lx_){ neighbourg(5) = i-Lx_; }
	else { neighbourg(5) = this->n_-Lx_+i; }

	return neighbourg;
}

template<typename Type>
double Triangle<Type>::occupation_number(Vector<double>& ni){
	double max(0);
	for(unsigned int i(0);i<Lx_;i++){
		for(unsigned int j(0);j<Ly_;j++){
			for(unsigned int k(0);k<this->m_;k++){
				ni(i+j*Lx_) += norm_squared(this->T_(i+j*Lx_,k));
			}
			if(ni(i+j*Lx_) > max){max = ni(i+j*Lx_);}
		}
	}
	return max;
}
#endif
