#include "SquareDimerizedBar.hpp"

SquareDimerizedBar::SquareDimerizedBar(System const& s, Vector<double> const& t):
	System(s),
	Square<double>(set_ab(),8,"square-dimerizedbar"),
	t_(t)
{
	if(t_.size()==8){
		if(status_==3){ init_lattice(); }
		if(status_==2){
			init_fermionic();
			same_wf_ = false;

			system_info_.text("SquareDimerizedBar :");
			system_info_.item("Each pair of color has a different Hamiltonian.");
			system_info_.item("8 sites in a 4x2 unit cell");
			system_info_.item("The dimerization occurs on bonds stacked in columns with strongest "+RST::math("\\vert t\\vert"));
			system_info_.item(RST::math("1=t_0\\geq t_i\\geq 0")+" (signs are set in H)");

			filename_ += "-t";
			for(unsigned int i(0);i<t_.size();i++){
				filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
			}
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : t must contain 8 values (currently contains "<<t_.size()<<")"<<std::endl; }
}

/*{method needed for running*/
void SquareDimerizedBar::compute_H(unsigned int const& c){
	H_.set(n_,n_,0);

	unsigned int t(0);
	int s(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		if(c<2){
			switch(obs_[0](i,5)){
				case 0: { t = obs_[0](i,3)?1:0; s =-(obs_[0](i,4)?bc_:1); }break;
				case 1: { t = obs_[0](i,3)?3:2; s = (obs_[0](i,4)?bc_:1); }break;
				case 2: { t = obs_[0](i,3)?1:4; s =-(obs_[0](i,4)?bc_:1); }break;
				case 3: { t = obs_[0](i,3)?3:2; s = (obs_[0](i,4)?bc_:1); }break;
				case 4: { t = obs_[0](i,3)?5:4; s = (obs_[0](i,4)?bc_:1); }break;
				case 5: { t = obs_[0](i,3)?6:7; s = (obs_[0](i,4)?bc_:1); }break;
				case 6: { t = obs_[0](i,3)?5:0; s = (obs_[0](i,4)?bc_:1); }break;
				case 7: { t = obs_[0](i,3)?6:7; s = (obs_[0](i,4)?bc_:1); }break;
			}
		} else {
			switch(obs_[0](i,5)){
				case 0: { t = obs_[0](i,3)?1:4; s =-(obs_[0](i,4)?bc_:1); }break;
				case 1: { t = obs_[0](i,3)?3:2; s = (obs_[0](i,4)?bc_:1); }break;
				case 2: { t = obs_[0](i,3)?1:0; s =-(obs_[0](i,4)?bc_:1); }break;
				case 3: { t = obs_[0](i,3)?3:2; s = (obs_[0](i,4)?bc_:1); }break;
				case 4: { t = obs_[0](i,3)?5:0; s = (obs_[0](i,4)?bc_:1); }break;
				case 5: { t = obs_[0](i,3)?6:7; s = (obs_[0](i,4)?bc_:1); }break;
				case 6: { t = obs_[0](i,3)?5:4; s = (obs_[0](i,4)?bc_:1); }break;
				case 7: { t = obs_[0](i,3)?6:7; s = (obs_[0](i,4)?bc_:1); }break;
			}
		}
		H_(obs_[0](i,0),obs_[0](i,1)) = s*t_(t);
	}
	H_ += H_.transpose();
}

void SquareDimerizedBar::create(){
	for(unsigned int c(0);c<N_;c++){
		if(!(c%2)){
			status_=2;
			compute_H(c);
			diagonalize(true);
		}
		if(status_==1){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
		}
	}
}

void SquareDimerizedBar::save_param(IOFiles& w) const {
	if(w.is_binary()){
		w<<t_;
		w.add_to_header()->title("t=("+my::tostring(t_)+")",'<');
		w.add_to_header()->add(system_info_.get());
	} else { w<<t_<<" "; }
}

Matrix<double> SquareDimerizedBar::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 4.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 0.0;
	tmp(1,1) = 2.0;
	return tmp;
}

unsigned int SquareDimerizedBar::unit_cell_index(Vector<double> const& x) const {
	if(my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.0, eq_prec_,eq_prec_)){ return 0; }
		if(my::are_equal(x(0),0.25,eq_prec_,eq_prec_)){ return 1; }
		if(my::are_equal(x(0),0.5, eq_prec_,eq_prec_)){ return 2; }
		if(my::are_equal(x(0),0.75,eq_prec_,eq_prec_)){ return 3; }
	} else {
		if(my::are_equal(x(0),0.0, eq_prec_,eq_prec_)){ return 4; }
		if(my::are_equal(x(0),0.25,eq_prec_,eq_prec_)){ return 5; }
		if(my::are_equal(x(0),0.5, eq_prec_,eq_prec_)){ return 6; }
		if(my::are_equal(x(0),0.75,eq_prec_,eq_prec_)){ return 7; }
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void SquareDimerizedBar::display_results(){
	compute_H(0);
	std::string t(my::tostring(t_));
	draw_lattice(true,true,false,(dir_nn_[2]+dir_nn_[3])*0.5,"-d:t "+t, "t=("+t+")");
}

void SquareDimerizedBar::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-dimerizedbar";
	display_results();

	//plot_band_structure();
}
/*}*/
