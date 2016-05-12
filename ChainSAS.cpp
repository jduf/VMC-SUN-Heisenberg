#include "ChainSAS.hpp"

ChainSAS::ChainSAS(System const& s, Vector<double> const& t):
	System(s),
	Chain<double>(set_spuc(N_),"chain-sas"),
	t_(t)
{
	std::cout<<J_<<std::endl;
	if(m_ == 1 && std::abs(J_(1)/J_(0))>0.2){ status_ = 4; }
	if(status_==2){
		init_fermionic();

		system_info_.text("ChainSAS:");
		system_info_.item("Work in the strong coupling limit");
		if(J_(0)<0.0){
			same_wf_ = false;
			system_info_.item("Each color has a different Hamiltonian.");
			system_info_.item("Feromagnetic strong coupling (symmetric irrep)");
		} else {
			same_wf_ = true;
			system_info_.item("Each color has the same Hamiltonian.");
			system_info_.item("Antiferomagnetic strong coupling (antisymmetric irrep)");
		}
		if(J_(1)<0.0){ system_info_.item("Feromagnetic week coupling "); }
		else         { system_info_.item("Antiferomagnetic week coupling"); }


		//filename_ += "-t"+std::string(t>=0?"+":"")+my::tostring(t_);
	}
}

/*{method needed for running*/
void ChainSAS::compute_H(unsigned int const& c){
	H_.set(n_,n_,0);

	if(c%2){
		for(unsigned int i(0);i<obs_[0].nlinks();i++){
			H_(obs_[0](i,0),obs_[0](i,1)) = (obs_[0](i,4)?bc_*t_(obs_[0](i,5)):t_(obs_[0](i,5)));
		}
	} else {
		for(unsigned int i(0);i<obs_[0].nlinks();i++){
			H_(obs_[0](i,0),obs_[0](i,1)) = -(obs_[0](i,4)?bc_*t_(obs_[0](i,5)):t_(obs_[0](i,5)));
		}
	}
	H_ += H_.transpose();
}

void ChainSAS::create(){
	if(same_wf_){
		compute_H(0);
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
	} else {
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
}

void ChainSAS::save_param(IOFiles& w) const {
	if(w.is_binary()){
		std::string s("t=("+my::tostring(t_)+")");

		w.add_header()->title(s,'<');
		w<<t_;
		w.add_header()->add(system_info_.get());
	} else { w<<t_<<" "; }
}

unsigned int ChainSAS::set_spuc(unsigned int const N){
	if(!(N%2)){ return N; }
	else {
		std::cerr<<__PRETTY_FUNCTION__<<" : works only for k=2"<<std::endl;
		return 0;
	}
}
/*}*/

/*{method needed for checking*/
void ChainSAS::check(){

}
/*}*/
