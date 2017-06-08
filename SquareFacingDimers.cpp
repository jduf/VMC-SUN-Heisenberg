#include "SquareFacingDimers.hpp"

SquareFacingDimers::SquareFacingDimers(System const& s, Vector<double> const& t):
	System(s),
	Square<double>(set_ab(),2,"square-facingdimers"),
	t_(t)
{
	if(t_.size()==3){
		if(status_==3){ init_lattice(); }
		if(status_==2){
			init_fermionic();

			system_info_.text("SquareFacingDimers :");
			system_info_.item("2 sites in a 2x1 unit cell");
			system_info_.item("Dimerization in columns");

			filename_ += "-t";
			for(unsigned int i(0);i<t_.size();i++){
				filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
			}
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : t must contain 3 values (currently contains "<<t_.size()<<")"<<std::endl; }
}

/*{method needed for running*/
void SquareFacingDimers::compute_H(){
	H_.set(n_,n_,0);

	double t(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		if(obs_[0](i,5)) { t = obs_[0](i,3)?t_(1):t_(2); }
		else{ t = obs_[0](i,3)?-t_(1):t_(0); }
		H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_:1)*t;
	}
	H_ += H_.transpose();
}

void SquareFacingDimers::create(){
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

void SquareFacingDimers::save_param(IOFiles& w) const {
	if(w.is_binary()){
		w<<t_;
		w.add_to_header()->title("t=("+my::tostring(t_)+")",'<');
		w.add_to_header()->add(system_info_.get());
	} else { w<<t_<<" "; }
}

Matrix<double> SquareFacingDimers::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 0.0;
	tmp(1,1) = 1.0;
	return tmp;
}

unsigned int SquareFacingDimers::unit_cell_index(Vector<double> const& x) const {
	if(my::are_equal(x(0),0.0,eq_prec_,eq_prec_)){ return 0; }
	if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_)){ return 1; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void SquareFacingDimers::display_results(){
	compute_H();
	std::string t(my::tostring(t_));
	draw_lattice(true,true,false,(dir_nn_[2]+dir_nn_[3])*0.5,"-d:t "+t, "t=("+t+")");
}

void SquareFacingDimers::check(){
	info_ = "";
	path_ = "";
	dir_ = "./";
	filename_ ="square-facingdimers";
	display_results();

	//plot_band_structure();
}

std::string SquareFacingDimers::extract_level_9(){
	Gnuplot gp(analyse_+path_+dir_,filename_);
	gp+="splot '"+filename_+".dat' u 2:3:4 w e notitle";
	gp.save_file();
	gp.create_image(true,"png");
	return filename_;
}
/*}*/
