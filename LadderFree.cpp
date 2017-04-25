#include "LadderFree.hpp"

LadderFree::LadderFree(System const& s, Vector<double> const& t, Vector<double> const& mu):
	System(s),
	Ladder<double>(set_spuc(t,mu,N_/m_),"ladder-free"),
	t_(t),
	mu_(mu)
{
	if(status_==2 && t_.ptr()){
		init_fermionic();

		system_info_.text("LadderFree :");
		system_info_.text(" Each color has the same Hamiltonian.");

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
void LadderFree::compute_H(){
	H_.set(n_,n_,0);

	unsigned int s0(0);
	unsigned int s1(0);
	unsigned int k(0);
	unsigned int l(0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		s0 = obs_[0](i,0);
		s1 = obs_[0](i,1);

		H_(s0,s1) = (obs_[0](i,4)?bc_*t_(k):t_(k));
		k = (k+1)%(3*spuc_/2);

		if(!obs_[0](i,3)){
			H_(s0,s0) = mu_(l)/2.0;
			l = (l+1)%spuc_;
		}
	}
	H_ += H_.transpose();
}

void LadderFree::create(){
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

void LadderFree::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("t=(");
		Vector<double> param(t_.size()+mu_.size());

		for(unsigned int i(0);i<t_.size()-1;i++){
			param(i) = t_(i);
			s += my::tostring(t_(i))+",";
		}
		param(t_.size()-1) = t_.back();
		s += my::tostring(t_.back())+") "+RST::math("\\mu")+"=(";

		for(unsigned int i(0);i<mu_.size()-1;i++){
			param(i+t_.size()) = mu_(i);
			s += my::tostring(mu_(i))+",";
		}
		param.back() = mu_.back();
		s += my::tostring(mu_.back())+")";

		w.add_to_header()->title(s,'<');
		w<<param;
		w.add_to_header()->add(system_info_.get());
	} else { w<<t_<<" "<<mu_<<" "; }
}

unsigned int LadderFree::set_spuc(Vector<double> const& t, Vector<double> const& mu, unsigned int const& spuc){
	if((t.size()*2/3)%spuc == 0 && mu.size()%spuc==0 && mu.size()<9){ return mu.size(); }
	else {
		std::cerr<<__PRETTY_FUNCTION__<<" : invalid or incoherent t and mu sizes : t:="<<t.size()<<", mu:="<<mu.size()<<std::endl;
		return 1;
	}
}
/*}*/

/*{method needed for checking*/
void LadderFree::display_results(){
	compute_H();
	std::string t(my::tostring(t_));
	std::string mu(my::tostring(mu_));
	draw_lattice(true,true,true,"ladder-free -d:t"+t+" -d:mu "+mu, "t=("+t+")"+RST::math("\\mu")+"=("+mu+")");
}

void LadderFree::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="ladder-free";
	display_results();

	//compute_H();
	//plot_band_structure();
	//create_obs(0);
	//std::cout<<obs_[3].get_links()<<std::endl;
}
/*}*/

/*{method needed for analysing*/
std::string LadderFree::extract_level_6(){
	(*data_write_)<<N_<<" "<<m_<<" "<<n_<<" "<<bc_<<" "<<asin(J_(1))<<" "<<obs_[0][0]<<IOFiles::endl;

	display_results();

	save_param(*jd_write_);
	save(*jd_write_);

	return filename_;
}
/*}*/
