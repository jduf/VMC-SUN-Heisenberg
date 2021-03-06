#include "Honeycomb0pp.hpp"

Honeycomb0pp::Honeycomb0pp(System const& s, double const& td):
	System(s),
	Honeycomb<double>(set_ab(),6,"honeycomb-0pp"),
	td_(td)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		system_info_.text("Honeycomb0pp :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("6 sites per unit cell");
		system_info_.item("1/3 of the hexagons with 0-flux.");
		system_info_.item("If td<0, each 0-flux hexagon is surrounded by pi-flux hexagons, 0-flux otherwise.");
		system_info_.item("th is set to -1");

		filename_ += "-td" + my::tostring(td_);
	}
}

/*{method needed for running*/
void Honeycomb0pp::compute_H(){
	H_.set(n_,n_,0);

	double th(-1.0);
	unsigned int s0(0);
	unsigned int s1(0);
	unsigned int ab0(0);
	unsigned int ab1(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		ab0 = obs_[0](i,5);
		ab1 = obs_[0](i,6);
		if((ab0==0 && ab1==1) || (ab0==2 && ab1==3) || (ab0==4 && ab1==5)){ H_(s0,s1) = (obs_[0](i,4)?bc_*td_:td_); }
		else { H_(s0,s1) = (obs_[0](i,4)?bc_*th:th); }
	}
	H_ += H_.transpose();
}

void Honeycomb0pp::create(){
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

void Honeycomb0pp::save_param(IOFiles& w) const {
	if(w.is_binary()){
		w<<Vector<double>(1,td_);
		w.add_to_header()->title(RST::math("t_d")+"="+my::tostring(td_)+", "+RST::math("t_h")+"=-1",'<');
		w.add_to_header()->add(system_info_.get());
	} else { w<<td_<<" "; }
}

Matrix<double> Honeycomb0pp::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 3.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 1.5;
	tmp(1,1) = 1.5*sqrt(3.0);
	return tmp;
}

unsigned int Honeycomb0pp::unit_cell_index(Vector<double> const& x) const {
	Vector<double> match(2,0);
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 0; }
	match(0) = 1.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 1; }
	match(1) = 1.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 2; }
	match(0) = 0;
	match(1) = 2.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 3; }
	match(0) = 2.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 4; }
	match(1) = 1.0/3.0;
	if(my::are_equal(x,match,eq_prec_,eq_prec_)){ return 5; }
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void Honeycomb0pp::display_results(){
	compute_H();
	std::string td(my::tostring(td_));
	draw_lattice(true,true,false,ref_(3)?(dir_nn_[4]+dir_nn_[3])*1.5:dir_nn_[3]*1.25+dir_nn_[4]*0.25,"-d:td "+td,RST::math("0\\pi\\pi")+", "+RST::math("t_d")+"="+td);
}

void Honeycomb0pp::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="honeycomb-0pp";
	display_results();

	//compute_H();
	//plot_band_structure();


	Data<double> b1;
	b1.merge(obs_[1][0]);
	b1.merge(obs_[1][4]);
	b1.merge(obs_[1][8]);
	b1.complete_analysis(1e-5);
	Data<double> b2;
	b2.merge(obs_[1][1]);
	b2.merge(obs_[1][2]);
	b2.merge(obs_[1][3]);
	b2.merge(obs_[1][5]);
	b2.merge(obs_[1][6]);
	b2.merge(obs_[1][7]);
	b2.complete_analysis(1e-5);

	std::cerr<<N_<<" "<<m_<<" "<<n_<<" "<<bc_<<" "<<td_<<" "<<obs_[0][0]<<" "<<b1<<" "<<b2<<" "<<ref_<<std::endl;
}
/*}*/
