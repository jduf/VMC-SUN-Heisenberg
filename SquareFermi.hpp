#ifndef DEF_SQUAREFERMI
#define DEF_SQUAREFERMI

#include "Square.hpp"

template<typename Type>
class SquareFermi: public Square<Type>{
	public:
		SquareFermi(System const& s);
		~SquareFermi() = default;

		void create();
		void check();

	protected:
		void compute_H();
		void display_results();

		/*!Set the unit cell's vectors*/
		Matrix<double> set_ab() const;
		/*!Returns the index of the site at position x in the unit cell*/
		unsigned int match_pos_in_ab(Vector<double> const& x) const { (void)(x); return 0; }
};

template<typename Type>
SquareFermi<Type>::SquareFermi(System const& s):
	System(s),
	Square<Type>(set_ab(),1,"square-fermi")
{
	if(this->status_==2){
		this->init_fermionic();

		this->system_info_.text("Fermi : all colors experience the same Hamiltonian");
	}
}

/*{method needed for running*/
template<typename Type>
void SquareFermi<Type>::compute_H(){
	this->H_.set(this->n_,this->n_,0);
	for(unsigned int i(0);i<this->obs_[0].nlinks(); i++){
		this->H_(this->obs_[0](i,0),this->obs_[0](i,1)) = (this->obs_[0](i,4)?this->bc_:1);
	}
	this->H_ += this->H_.transpose();
}

template<typename Type>
Matrix<double> SquareFermi<Type>::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1;
	tmp(1,0) = 0;
	tmp(0,1) = 0;
	tmp(1,1) = 1;
	return tmp;
}
/*}*/

/*{method needed for checking*/
template<typename Type>
void SquareFermi<Type>::display_results(){
	compute_H();

	std::string color("black");
	std::string linestyle("solid");
	Vector<double> xy0(2,0);
	Vector<double> xy1(2,0);
	PSTricks ps(this->info_+this->path_+this->dir_,this->filename_);
	ps.begin(-20,-20,20,20,this->filename_);
	ps.polygon(this->lattice_corners_,"linecolor=green");

	double x_shift(-this->ab_(0,0)/2.0);
	double y_shift(-this->ab_(1,1)/2.0);
	Matrix<double> polygon(4,2);
	polygon(0,0)=x_shift;
	polygon(0,1)=y_shift;
	polygon(1,0)=x_shift+this->ab_(0,0);
	polygon(1,1)=y_shift+this->ab_(1,0);
	polygon(2,0)=x_shift+this->ab_(0,0)+this->ab_(0,1);
	polygon(2,1)=y_shift+this->ab_(1,0)+this->ab_(1,1);
	polygon(3,0)=x_shift+this->ab_(0,1);
	polygon(3,1)=y_shift+this->ab_(1,1);
	ps.polygon(polygon,"linecolor=black");

	Type t;
	unsigned int s0;
	unsigned int s1;
	for(unsigned int i(0);i<this->obs_[0].nlinks();i++) {
		s0 = this->obs_[0](i,0);
		xy0 = this->x_[s0];

		s1 = this->obs_[0](i,1);
		xy1 = this->x_[s1];

		t = this->H_(s0,s1);
		if(std::abs(t)>1e-4){
			if((xy0-xy1).norm_squared()>1.0001){
				linestyle = "dashed";
				xy1 = (xy0+this->dir_nn_[this->obs_[0](i,3)]).chop();
				ps.put(xy1(0)-0.20,xy1(1)+0.15,"\\tiny{"+my::tostring(s1)+"}");
			} else{ linestyle = "solid"; }

			if(my::real(t)<0){ color = "red"; }
			else { color = "blue"; }
			ps.line("-",xy0(0),xy0(1),xy1(0),xy1(1), "linewidth=1pt,linecolor="+color+",linestyle="+linestyle);
		}
		if(i%2){ ps.put(xy0(0)-0.20,xy0(1)+0.15,"\\tiny{"+my::tostring(s0)+"}"); }
	}
	ps.line("-",this->boundary_[0](0),this->boundary_[0](1),this->boundary_[1](0),this->boundary_[1](1),"linecolor=yellow");
	ps.line("-",this->boundary_[3](0),this->boundary_[3](1),this->boundary_[0](0),this->boundary_[0](1),"linecolor=yellow");
	ps.end(true,true,true);
}

template<typename Type>
void SquareFermi<Type>::check(){
	this->info_ = "";
	this->path_ = "";
	this->dir_  = "./";
	this->filename_ ="square-fermi";
	this->display_results();
}
/*}*/
#endif
