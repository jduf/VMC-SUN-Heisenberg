#include "SquareFree.hpp"

SquareFree::SquareFree(System const& s, Vector<double> const& t, Vector<double> const& mu):
	System(s),
	Square<double>(set_ab(ref_(3)),6,"square-free"),
	t_(t),
	mu_(mu)
{
	if(status_==3){ init_lattice(); }
	if(status_==2){
		init_fermionic();

		//bool need_compute_additional_links(true);
		//for(unsigned int i(0);i<obs_.size();i++){
			//if(obs_[i].get_type() == 4){ need_compute_additional_links = false; i=obs_.size(); }
		//}
		//if(need_compute_additional_links){ init_additional_links(); }

		same_wf_ = false;

		system_info_.text("SquareFree :");
		system_info_.item("Each color has a different Hamiltonian.");
		//system_info_.item("There is an additional second neighbour hopping for every 1/k sites.");

		filename_ += "-t";
		for(unsigned int i(0);i<t_.size();i++){
			filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
		}
		filename_ += "-mu";
		for(unsigned int i(0);i<mu_.size();i++){
			filename_ += ((mu_(i)>=0)?"+":"")+my::tostring(mu_(i));
		}
	}
}

/*{method needed for running*/
void SquareFree::init_additional_links(){
	Vector<double> x;
	Matrix<int> tmp(n_*2,3);
	for(unsigned int i(0);i<n_;i++){
		x = x_[i]+dir_nn_[0]+dir_nn_[1]*2.0;
		tmp(2*i,2) = cross_boundary(x_[i],x);
		tmp(2*i,0) = i;
		tmp(2*i,1) = site_index(x);

		x = x_[i]-dir_nn_[1]+dir_nn_[0]*2.0;
		tmp(2*i+1,2) = cross_boundary(x_[i],x);
		tmp(2*i+1,0) = i;
		tmp(2*i+1,1) = site_index(x);
	}
	obs_.push_back(Observable("Additional links",4,0,tmp));
}

void SquareFree::compute_H(unsigned int const& c){
	H_.set(n_,n_,0);

	unsigned int s0(0);
	unsigned int s1(0);
	unsigned int ab(0);
	unsigned int l(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);
		ab = obs_[0](i,5);
		l = (2*ab+obs_[0](i,3));
		if(c<2){
			switch(l){
				case 0:{  H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(4); }break;
				case 4:{  H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(0); }break;
				case 8:{  H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(12); }break;
				case 12:{ H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(8); }break;
				default:{ H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(l); }
			}
		} else {
			switch(l){
				case 0:{  H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(0); }break;
				case 4:{  H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(4); }break;
				case 8:{  H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(8); }break;
				case 12:{ H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(12); }break;
				default:{ H_(s0,s1) = (obs_[0](i,4)?bc_:1)*t_(l); }
			}
		}
		if(obs_[0](i,3)){ H_(s0,s0) = mu_((ab+c)%mu_.size())/2; }
	}
	H_ += H_.transpose();
}

void SquareFree::create(){
	for(unsigned int c(0);c<N_;c++){
		status_ = 2;
		compute_H(c);
		diagonalize(true);
		if(status_==1){
			for(unsigned int i(0);i<n_;i++){
				for(unsigned int j(0);j<M_(c);j++){
					EVec_[c](i,j) = H_(i,j);
				}
			}
		}
	}
}

void SquareFree::save_param(IOFiles& w) const {
	if(w.is_binary()){
		Vector<double> param(t_.size()+mu_.size());
		for(unsigned int i(0);i<t_.size();i++){ param(i) = t_(i); }
		for(unsigned int i(0);i<mu_.size();i++){ param(i+t_.size()) = mu_(i); }
		w<<param;
		w.add_to_header()->title(RST::math("\\mu")+"=("+my::tostring(mu_)+") t=("+my::tostring(t_)+")",'<');
		w.add_to_header()->add(system_info_.get());
	} else { w<<t_<<" "<<mu_<<" "; }
}

Matrix<double> SquareFree::set_ab(unsigned int const& ref3) const {
	Matrix<double> tmp(2,2);
	//if(ref3==2){
	//tmp(0,0) = 2.0;
	//tmp(1,0) = 1.0;
	//tmp(0,1) =-1.0;
	//tmp(1,1) = 2.0;
	//} else {
	//tmp(0,0) = 2.0;
	//tmp(1,0) =-1.0;
	//tmp(0,1) = 1.0;
	//tmp(1,1) = 2.0;
	//}
	//(void)(ref3);
	//tmp(0,0) = 2.0;
	//tmp(1,0) = 0.0;
	//tmp(0,1) = 0.0;
	//tmp(1,1) = 2.0;
	
	(void)(ref3);
	tmp(0,0) = 4.0;
	tmp(1,0) = 0.0;
	tmp(0,1) = 0.0;
	tmp(1,1) = 2.0;
	
	return tmp;
}

unsigned int SquareFree::unit_cell_index(Vector<double> const& x) const {
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
void SquareFree::display_results(){
	compute_H(0);
	std::string t(my::tostring(t_));
	std::string mu(my::tostring(mu_));
	draw_lattice(true,true,false,(dir_nn_[2]+dir_nn_[3])*0.5,"-d:t "+t+" -d:mu "+mu,"t=("+t+") "+RST::math("\\mu")+"=("+mu+")");
}

void SquareFree::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="square-free";
	display_results();

	//compute_H();
	//plot_band_structure();
}
/*}*/
