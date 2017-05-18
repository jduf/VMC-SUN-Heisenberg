#include "LadderDimerA.hpp"

LadderDimerA::LadderDimerA(System const& s, Vector<double> const& t):
	System(s),
	Ladder<double>(2,"ladder-dimerA"),
	t_(t)
{
	if(t_.size()==2){
		if(status_==2){
			init_fermionic();

			system_info_.text("LadderDimerA :");
			system_info_.text(" Each color has the same Hamiltonian.");
			system_info_.text(" Dimer");

			filename_ += "-t";
			for(unsigned int i(0);i<t_.size();i++){
				filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
			}
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : t must contain 2 values (currently contains "<<t_.size()<<")"<<std::endl; }
}

/*{method needed for running*/
void LadderDimerA::compute_H(){
	H_.set(n_,n_,0);

	double t(0.0);
	for(unsigned int i(0);i<obs_[0].nlinks();i++){
		if(obs_[0](i,5)){ t = t_(0); }
		else {
			if(obs_[0](i,3)){ t = t_(1); }
			else            { t = t_(0); }
		}

		H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t:t);
	}
	H_ += H_.transpose();
}

void LadderDimerA::create(){
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

void LadderDimerA::save_param(IOFiles& w) const {
	if(w.is_binary()){
		w<<t_;
		w.add_to_header()->title("t=("+my::tostring(t_)+")",'<');
		w.add_to_header()->add(system_info_.get());
	} else { w<<t_<<" "; }
}
/*}*/

/*{method needed for checking*/
void LadderDimerA::display_results(){
	compute_H();
	std::string t(my::tostring(t_));
	draw_lattice(true,true,true,"-d:t "+t, "t=("+t+")");
}

void LadderDimerA::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="ladder-dimerA";
	display_results();

	//compute_H();
	//plot_band_structure();
}
/*}*/
