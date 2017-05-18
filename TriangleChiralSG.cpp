#include "TriangleChiralSG.hpp"

TriangleChiralSG::TriangleChiralSG(System const& s, double const& phi):
	System(s),
	Triangle<std::complex<double> >(set_ab(),1,"triangle-chiralSG"),
	phi_(phi)
{
	if(1.0*phi_<=N_/m_){
		if(status_==3){ init_lattice(); }
		if(status_==2){
			init_fermionic();

			system_info_.text("TriangleChiralSG :");
			system_info_.item("Each color has the same Hamiltonian.");
			system_info_.item("There is a flux of "+RST::math(my::tostring(phi)+"\\pi/"+my::tostring(N_/m_)) + "per triangular plaquette.");
			system_info_.item("String gauge (no proper unit cell)");

			filename_ += "-phi"+my::tostring(phi_);
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : the flux per plaquette shouldn't be bigger than pi"<<std::endl; }
}

/*{method needed for running*/
void TriangleChiralSG::compute_H(){
	Matrix<unsigned int> nn(n_,6);
	for(unsigned int i(0);i<n_;i++){
		Matrix<int> tmp(get_neighbourg(i));
		for(unsigned int n(0);n<6;n++){ nn(i,n) = tmp(n,0); }
	}

	Matrix<double> M(2,2);
	M(0,0) = dir_nn_[0](0);
	M(1,0) = dir_nn_[0](1);
	M(0,1) = dir_nn_[1](0);
	M(1,1) = dir_nn_[1](1);
	Lapack<double>(M,false,'G').inv();

	H_.set(n_,n_,0);
	Vector<double> move;
	unsigned int s;
	unsigned int m;
	for(unsigned int i(0);i<n_;i++){
		move = (M*(x_[i]-x_[0])).chop();
		s = 0;
		if(move(0)>0){
			m = std::round(move(0));
			for(unsigned int j(0);j<m;j++){
				H_(s,nn(s,1)) -= 2.0;
				H_(nn(s,1),s) += 2.0;
				s = nn(s,0);
				H_(s,nn(s,2)) -= 2.0;
				H_(nn(s,2),s) += 2.0;
			}
		}
		if(move(0)<0){
			m = std::round(-move(0));
			for(unsigned int j(0);j<m;j++){
				H_(s,nn(s,2)) += 2.0;
				H_(nn(s,2),s) -= 2.0;
				s = nn(s,3);
				H_(s,nn(s,1)) += 2.0;
				H_(nn(s,1),s) -= 2.0;
			}
		}
		if(move(1)>0){
			m = std::round(move(1));
			for(unsigned int j(0);j<m;j++){
				s = nn(s,1);
				H_(s,nn(s,2)) -= 2.0;
				H_(s,nn(s,3)) -= 2.0;
				H_(nn(s,2),s) += 2.0;
				H_(nn(s,3),s) += 2.0;
			}
		}
		if(move(1)<0){
			m = std::round(-move(1));
			for(unsigned int j(0);j<m;j++){
				H_(s,nn(s,2)) += 2.0;
				H_(s,nn(s,3)) += 2.0;
				H_(nn(s,2),s) -= 2.0;
				H_(nn(s,3),s) -= 2.0;
				s = nn(s,4);
			}
		}
		H_(s,nn(s,1)) -= 1.0;
		H_(nn(s,1),s) += 1.0;
	}

	int p;
	double t(-1.0);
	double phi(phi_*m_*M_PI/N_);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		p = std::real(H_(obs_[0](i,0),obs_[0](i,1)));
		p = p%int(2*N_/m_);
		H_(obs_[0](i,0),obs_[0](i,1)) = std::polar((obs_[0](i,4)?t*bc_:t),p*phi);

		p = std::real(H_(obs_[0](i,1),obs_[0](i,0)));
		p = p%int(2*N_/m_);
		H_(obs_[0](i,1),obs_[0](i,0)) = std::polar((obs_[0](i,4)?t*bc_:t),p*phi);
	}

}

void TriangleChiralSG::create(){
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

void TriangleChiralSG::save_param(IOFiles& w) const {
	if(w.is_binary()){
		w<<Vector<double>(1,phi_);
		w.add_to_header()->title(RST::math("\\phi")+"="+my::tostring(phi_),'<');
		w.add_to_header()->add(system_info_.get());
	} else { w<<phi_<<" "; }
}

Matrix<double> TriangleChiralSG::set_ab() const {
	Matrix<double> tmp(2,2);
	tmp(0,0) = 1;
	tmp(1,0) = 0.0;
	tmp(0,1) = 0.5;
	tmp(1,1) = sqrt(3.0)/2.0;
	return tmp;
}

void TriangleChiralSG::bond_energy_obs(){
	unsigned int idx(obs_.size());
	obs_.push_back(Observable("Bond energy",1,obs_[0].nlinks(),obs_[0].nlinks()));
	obs_[idx].remove_links();

	for(unsigned int i(0);i<obs_[0].nlinks();i++){ obs_[0](i,2) = i; }
}
/*}*/

/*{method needed for checking*/
void TriangleChiralSG::display_results(){
	compute_H();
	std::string phi(my::tostring(phi_));
	draw_lattice(false,true,false,dir_nn_[3]*1.75+dir_nn_[4]*0.25,"-d:phi "+phi,RST::math("\\phi")+"="+phi);
}

void TriangleChiralSG::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="triangle-chiralSG";
	display_results();

	//compute_H();
	//std::cout<<H_<<std::endl;

	//compute_H();
	//plot_band_structure();
}
/*}*/
