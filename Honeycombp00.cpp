#include "Honeycombp00.hpp"

Honeycombp00::Honeycombp00(System const& s, double const& td):
	System(s),
	Honeycomb<double>(set_ab(),12,"honeycomb-p00"),
	td_(td)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("Honeycombp00 :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("6 sites per unit cell");
		system_info_.item("1/3 of the hexagons with pi-flux.");
		system_info_.item("If td<0, each pi-flux hexagon is surrounded by 0-flux hexagons, pi-flux otherwise.");
		system_info_.item("th is set to -1");

		filename_ += "-td" + my::tostring(td_);
	}
}

/*{method needed for running*/
void Honeycombp00::compute_H(){
	H_.set(n_,n_,0);

	double th(-1.0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		switch(obs_[0](i,5)){
			case 0:
				{
					switch(obs_[0](i,6)){
						case 1: { H_(obs_[0](i,0),obs_[0](i,1)) =-(obs_[0](i,4)?bc_*th:th); }break;
						case 3: { H_(obs_[0](i,0),obs_[0](i,1)) =-(obs_[0](i,4)?bc_*th:th); }break;
						case 11:{ H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*td_:td_); }break;
					}
				}break;
			case 2:
				{
					switch(obs_[0](i,6)){
						case 3: { H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*td_:td_); }break;
						case 5: { H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*th:th); }break;
						case 1: { H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*th:th); }break;
					}
				}break;
			case 4:
				{
					switch(obs_[0](i,6)){
						case 5: { H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*th:th); }break;
						case 7: { H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*td_:td_); }break;
						case 3: { H_(obs_[0](i,0),obs_[0](i,1)) =-(obs_[0](i,4)?bc_*th:th); }break;
					}
				}break;
			case 6:
				{
					switch(obs_[0](i,6)){
						case 7: { H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*th:th); }break;
						case 9: { H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*th:th); }break;
						case 5: { H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*td_:td_); }break;
					}
				}break;
			case 8:
				{
					switch(obs_[0](i,6)){
						case 9: { H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*td_:td_); }break;
						case 11:{ H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*th:th); }break;
						case 7: { H_(obs_[0](i,0),obs_[0](i,1)) =-(obs_[0](i,4)?bc_*th:th); }break;
					}
				}break;
			case 10:
				{
					switch(obs_[0](i,6)){
						case 11:{ H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*th:th); }break;
						case 1: { H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*td_:td_); }break;
						case 9: { H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*th:th); }break;
					}
				}break;
		}
	}
	H_ += H_.transpose();
}

void Honeycombp00::create(){
	compute_H();
	diagonalize(true);
	for(unsigned int c(0);c<N_;c++){
		for(unsigned int i(0);i<n_;i++){
			for(unsigned int j(0);j<M_(c);j++){
				EVec_[c](i,j) = H_(i,j);
			}
		}
	}
}

void Honeycombp00::save_param(IOFiles& w) const {
	if(w.is_binary()){
		w<<Vector<double>(1,td_);
		w.add_to_header()->title(RST::math("t_d")+"="+my::tostring(td_)+", "+RST::math("t_h")+"=-1",'<');
		w.add_to_header()->add(system_info_.get());
	} else { w<<td_<<" "; }
}

Matrix<double> Honeycombp00::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 3.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 0.0;
	tmp(1,1) = 3.0*sqrt(3.0);
	return tmp;
}

unsigned int Honeycombp00::unit_cell_index(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x(0),0.0,eq_prec_,eq_prec_)){
		if(my::are_equal(x(1),0.0,    eq_prec_,eq_prec_)){ return 0; }
		if(my::are_equal(x(1),1.0/3.0,eq_prec_,eq_prec_)){ return 4; }
		if(my::are_equal(x(1),2.0/3.0,eq_prec_,eq_prec_)){ return 8; }
	}
	if(my::are_equal(x(0),1.0/3.0,eq_prec_,eq_prec_)){
		if(my::are_equal(x(1),0.0,    eq_prec_,eq_prec_)){ return 1; }
		if(my::are_equal(x(1),1.0/3.0,eq_prec_,eq_prec_)){ return 5; }
		if(my::are_equal(x(1),2.0/3.0,eq_prec_,eq_prec_)){ return 9; }
	}
	if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_)){
		if(my::are_equal(x(1),1.0/6.0,eq_prec_,eq_prec_)){ return 2; }
		if(my::are_equal(x(1),3.0/6.0,eq_prec_,eq_prec_)){ return 6; }
		if(my::are_equal(x(1),5.0/6.0,eq_prec_,eq_prec_)){ return 10;}
	}
	if(my::are_equal(x(0),5.0/6.0,eq_prec_,eq_prec_)){
		if(my::are_equal(x(1),1.0/6.0,eq_prec_,eq_prec_)){ return 3; }
		if(my::are_equal(x(1),3.0/6.0,eq_prec_,eq_prec_)){ return 7; }
		if(my::are_equal(x(1),5.0/6.0,eq_prec_,eq_prec_)){ return 11;}
	}

	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void Honeycombp00::display_results(){
	compute_H();
	std::string td(my::tostring(td_));
	draw_lattice(false,true,false,ref_(3)?(dir_nn_[3]+dir_nn_[4]+dir_nn_[5])*1.5:dir_nn_[3]*1.5+dir_nn_[4]+dir_nn_[5],"-d:td "+td,RST::math("\\pi 00")+" with "+RST::math("t_d")+"="+td);
}

void Honeycombp00::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="honeycomb-p00";
	//display_results();

	compute_H();
	plot_band_structure();
}
/*}*/
