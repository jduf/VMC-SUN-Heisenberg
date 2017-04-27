#include "SquareFermi.hpp"

SquareFermi::SquareFermi(System const& s):
	System(s),
	Square<double>(set_ab(),1,"square-fermi")
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("SquareFermi :");
		system_info_.item("Each color has the same Hamiltonian.");
	}
}

/*{method needed for running*/
void SquareFermi::compute_H(){
	H_.set(n_,n_,0);

	double t(-1.0);
	for(unsigned int i(0);i<obs_[0].nlinks(); i++){
		H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t:t);
	}
	H_ += H_.transpose();
}

void SquareFermi::create(){
	compute_H();
	diagonalize(true);
	if(status_==1){
		for(unsigned int c(0);c<N_;c++){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
		}
	}
}

Matrix<double> SquareFermi::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1;
	tmp(1,0) = 0;
	tmp(0,1) = 0;
	tmp(1,1) = 1;
	return tmp;
}
/*}*/

/*{method needed for checking*/
void SquareFermi::display_results(){
	compute_H();
	draw_lattice(false,true,false,(dir_nn_[2]+dir_nn_[3])*0.5,"","Fermi");
}

void SquareFermi::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-fermi";
	display_results();
}

void SquareFermi::param_fit_therm_limit(std::string& f, std::string& param, std::string& range){
	f="f(x) = a+b*x*x";
	param = "a,b";
	range = "[0:0.015]";
}
/*}*/
