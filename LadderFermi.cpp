#include"LadderFermi.hpp"

LadderFermi::LadderFermi(System const& s):
	System(s),
	Ladder(1,"ladder-fermi")
{
	if(status_==2){
		init_fermionic();

		system_info_.text("LadderFermi :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("Uniform real hopping term.");
	}
}

/*{method needed for running*/
void LadderFermi::compute_H(){
	H_.set(n_,n_,0);
	for(unsigned int i(0);i<obs_[0].nlinks(); i++){
		H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_:1);
	}
	H_ += H_.transpose();
}

void LadderFermi::create(){
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
/*}*/

/*{method needed for checking*/
void LadderFermi::display_results(){
	compute_H();
	draw_lattice(true,true,true,"","Fermi");
}

void LadderFermi::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="ladder-fermi";
	display_results();

	//plot_band_structure();
}
/*}*/
