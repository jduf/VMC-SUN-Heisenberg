#include "LadderRectangularPlaquetteA.hpp"

LadderRectangularPlaquetteA::LadderRectangularPlaquetteA(System const& s, Vector<double> const& t):
	System(s),
	Ladder<double>(6,"ladder-rectangularplaquetteA"),
	t_(t)
{
	if(t_.size()==3){
		std::cerr<<__PRETTY_FUNCTION__<<" : deprecated, t should contain 4 values"<<std::endl;
		if(status_==2){
			init_fermionic();

			system_info_.text("LadderRectangularPlaquetteA :");
			system_info_.text(" Each color has the same Hamiltonian.");
			system_info_.text(" Rectangular plaquette in a 6 site unit cell");
			system_info_.text(" pi flux in each square plaquette");

			filename_ += "-t";
			for(unsigned int i(0);i<t_.size();i++){
				filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
			}
		}
	} else if(t_.size()==4){
		if(status_==2){
			init_fermionic();

			system_info_.text("LadderRectangularPlaquetteA :");
			system_info_.text(" Each color has the same Hamiltonian.");
			system_info_.text(" Rectangular plaquette in a 6 site unit cell");
			system_info_.text(" pi flux in each square plaquette");

			filename_ += "-t";
			for(unsigned int i(0);i<t_.size();i++){
				filename_ += ((t_(i)>=0)?"+":"")+my::tostring(t_(i));
			}
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : t must contain 4 values (currently contains "<<t_.size()<<")"<<std::endl; }
}

/*{method needed for running*/
void LadderRectangularPlaquetteA::compute_H(){
	H_.set(n_,n_,0);

	double t(0.0);
	if(t_.size() == 3){// to preserve backward compatibility
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
						if(obs_[0](i,3)){ t = t_(1); }
						else            { t = t_(0); }
					}break;
				case 3:
					{ t = -t_(0); }break;
				case 4:
					{
						if(obs_[0](i,3)){ t = t_(1); }
						else            { t = t_(2); }
					}break;
				case 5:
					{ t = -t_(2); }break;
			}

			H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t:t);
		}
	} else {
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
				case 4:
					{
						if(obs_[0](i,3)){ t = t_(1); }
						else            { t = t_(3); }
					}break;
				case 5:
					{ t = -t_(3); }break;
			}

			H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t:t);
		}
	}
	H_ += H_.transpose();
}

void LadderRectangularPlaquetteA::create(){
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

void LadderRectangularPlaquetteA::save_param(IOFiles& w) const {
	if(w.is_binary()){
		w<<t_;
		w.add_to_header()->title("t=("+my::tostring(t_)+")",'<');
		w.add_to_header()->add(system_info_.get());
	} else { w<<t_<<" "; }
}
/*}*/

/*{method needed for checking*/
void LadderRectangularPlaquetteA::display_results(){
	compute_H();
	std::string t(my::tostring(t_));
	draw_lattice(true,true,true,"-d:t "+t, "t=("+t+")");
}

void LadderRectangularPlaquetteA::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="ladder-rectangularplaquetteA";
	display_results();

	//compute_H();
	//plot_band_structure();
}
/*}*/
