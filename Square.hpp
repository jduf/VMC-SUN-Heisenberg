#ifndef DEF_SQUARE
#define DEF_SQUARE

#include "CreateSystem.hpp"

template<typename Type>
class Square: public CreateSystem<Type>{
	public:
		Square(Parseur& P);
		~Square();

	protected:
		unsigned int Lx_;//!< dimension of the lattice along x-axis
		unsigned int Ly_;//!< dimension of the lattice along y-axis
		Matrix<Type> Px_;//!< translation operator along x-axis 
		Matrix<Type> Py_;//!< translation operator along y-axis 

		void compute_H();
		void study_system();
		void save(std::string filename);
};

template<typename Type>
Square<Type>::Square(Parseur& P):
	CreateSystem<Type>(P,4),
	Lx_(floor(sqrt(this->n_))),
	Ly_(floor(sqrt(this->n_)))
{
	this->bc_= P.get<double>("bc");
	if(!P.status()){
		if(this->n_==Ly_*Lx_){
			compute_H();
			this->compute_sts();
		} else {
			std::cerr<<"Square : the cluster is not a square"<<std::endl;
		}
	}
}

template<typename Type>
Square<Type>::~Square(){}

template<typename Type>
void Square<Type>::compute_H(){
	for(unsigned int i(0); i < this->n_; i++){
		/*horizontal hopping*/
		if( (i+1) % Lx_ ){ this->H_(i,i+1) = 1;}	
		else { this->H_(i+1-Lx_,i) = 1;}
		/*vertical hopping*/
		if( i+Lx_ < this->n_ ){  this->H_(i,i+Lx_) = 1; } 
		else { this->H_(i-(Ly_-1)*Lx_,i) = 1;}
	}
	this->H_ += this->H_.transpose();
}

template<typename Type>
void Square<Type>::study_system(){
	Matrix<double> ad_ia_i(Lx_,Ly_,0);

	for(unsigned int i(0);i<Ly_;i++){
		for(unsigned int j(0);j<Lx_;j++){
			for(unsigned int k(0);k<this->m_;k++){
				ad_ia_i(i,j) += norm_squared(this->T_(i+j*Lx_,k));
			}
		}
	}
	double max(0.);
	double min(0.);
	for(unsigned int i(0);i<Ly_;i++){
		for(unsigned int j(0);j<Lx_;j++){
			if(ad_ia_i(i,j) > max){ max = ad_ia_i(i,j); }
			if(ad_ia_i(i,j) < min){ min = ad_ia_i(i,j); }
		}
	}
	ad_ia_i -= min;
	ad_ia_i /= max;

	Vector<double> xy(Lx_);
	for(unsigned int i(0);i<Lx_;i++){
		xy(i) = i;
	}
	Gnuplot gp("AA","plot");
	gp.preplot("set cbrange [0:1]");
	gp.preplot("unset border");
	gp.preplot("unset xtics");
	gp.preplot("unset ytics");
	gp.preplot("unset key");
	gp.preplot("set size square");
	gp.save_data("surf",xy,xy,ad_ia_i);
	gp.add_plot_param("w image\n");
}
#endif

