#include "KagomePlaquette6A.hpp"

KagomePlaquette6A::KagomePlaquette6A(System const& s, double const& td):
	System(s),
	Kagome<double>(set_ab(),6,"kagome-plaquette6A"),
	td_(td)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		filename_ += ((td_>=0)?"-td+":"-td") + my::tostring(td_);

		system_info_.text("KagomePlaquette6A :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("6 sites per unit cell.");
		if(td_<0.0){ system_info_.item(RST::math("(0;0,\\pi)")+"-flux"); }
		if(td_>0.0){ system_info_.item(RST::math("(\\pi;\\pi,\\pi)")+"-flux"); }
	}
}

/*{method needed for running*/
void KagomePlaquette6A::compute_H(){
	H_.set(n_,n_,0);

	double th(-1.0);
	double t(0.0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		t=0.0;
		switch(obs_[0](i,5)){
			case 0: { t = th; } break;
			case 1: { t = td_; } break;
			case 2: { t = (obs_[0](i,3)?-th:td_); } break;
			case 3: { t = (obs_[0](i,3)?th:-th);  } break;
			case 4: { t = td_; } break;
			case 5: { t = (obs_[0](i,3)?th:td_); } break;
		}
		H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t:t);
	}
	H_ += H_.transpose();
}

void KagomePlaquette6A::create(){
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

void KagomePlaquette6A::save_param(IOFiles& w) const {
	if(w.is_binary()){
		w<<Vector<double>(1,td_);
		w.add_to_header()->title(RST::math("t_d")+"="+my::tostring(td_)+", "+RST::math("t_h")+"=-1",'<');
		w.add_to_header()->add(system_info_.get());
	} else { w<<td_<<" "; }
}

Matrix<double> KagomePlaquette6A::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 2.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 0.0;
	tmp(1,1) = 2.0*sqrt(3.0);
	return tmp;
}

unsigned int KagomePlaquette6A::unit_cell_index(Vector<double> const& x) const {
	if(my::are_equal(x(1),0.0,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.0,eq_prec_,eq_prec_)){ return 0; }
		if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_)){ return 2; }
	}
	if(my::are_equal(x(1),0.25,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.25,eq_prec_,eq_prec_)){ return 1; }
	}
	if(my::are_equal(x(1),0.5,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.0,eq_prec_,eq_prec_)){ return 5; }
		if(my::are_equal(x(0),0.5,eq_prec_,eq_prec_)){ return 3; }
	}
	if(my::are_equal(x(1),0.75,eq_prec_,eq_prec_)){
		if(my::are_equal(x(0),0.75,eq_prec_,eq_prec_)){ return 4; }
	}
	std::cerr<<__PRETTY_FUNCTION__<<" : unknown position in ab for x="<<x<<std::endl;
	return spuc_;
}
/*}*/

/*{method needed for checking*/
void KagomePlaquette6A::display_results(){
	compute_H();
	std::string td(my::tostring(td_));
	draw_lattice(false,true,false,(dir_nn_[2]+dir_nn_[3])*0.5,"-d:td "+td,RST::math("t_d")+"="+td);
}

void KagomePlaquette6A::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="kagome-plaquette6A";
	display_results();

	//compute_H();
	//plot_band_structure();
}
/*}*/
