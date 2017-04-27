#include "LadderSquarePlaquetteC.hpp"

LadderSquarePlaquetteC::LadderSquarePlaquetteC(System const& s, Vector<double> const& t):
	System(s),
	Ladder<double>(4,"ladder-squareplaquetteC"),
	t_(t)
{
	if(t_.size()==3){
		if(status_==2){
			init_fermionic();

			system_info_.text("LadderSquarePlaquetteC :");
			system_info_.text(" Each color has the same Hamiltonian.");
			system_info_.text(" Square plaquette in a 4 site unit cell");
			system_info_.text(" pi flux between the plaquettes");
			system_info_.text(" pi flux inside the plaquettes");

			filename_ += "-t";
			for(unsigned int i(0);i<t_.size();i++){
				filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
			}
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : t must contain 3 values (currently contains "<<t_.size()<<")"<<std::endl; }
}

/*{method needed for running*/
void LadderSquarePlaquetteC::compute_H(){
	H_.set(n_,n_,0);

	double t(0.0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		switch(obs_[0](i,5)){
			case 0:
				{
					if(obs_[0](i,3)){ t = t_(1); }
					else            { t = t_(0); }
				}break;
			case 1:
				{ t = -t_(0); }break;
			case 2:
				{
					if(obs_[0](i,3)){ t = t_(2); }
					else            { t = t_(0); }
				}break;
			case 3:
				{ t = -t_(0); }break;
		}

		H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t:t);
	}
	H_ += H_.transpose();
}

void LadderSquarePlaquetteC::create(){
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

void LadderSquarePlaquetteC::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("t=(");
		Vector<double> param(t_.size());

		for(unsigned int i(0);i<t_.size()-1;i++){
			param(i) = t_(i);
			s += my::tostring(t_(i))+",";
		}
		param(t_.size()-1) = t_.back();
		s += my::tostring(t_.back())+")";

		w.add_to_header()->title(s,'<');
		w<<param;
		w.add_to_header()->add(system_info_.get());
	} else { w<<t_<<" "; }
}
/*}*/

/*{method needed for checking*/
void LadderSquarePlaquetteC::display_results(){
	compute_H();
	std::string t(my::tostring(t_));
	draw_lattice(true,true,true,"-d:t "+t, "t=("+t+")");
}

void LadderSquarePlaquetteC::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="ladder-squareplaquetteC";
	display_results();

	//compute_H();
	//plot_band_structure();
}
/*}*/
