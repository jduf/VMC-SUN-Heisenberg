#include"LadderFreeFlux.hpp"

LadderFreeFlux::LadderFreeFlux(System const& s, Vector<double> const& t, Vector<double> const& flux):
	System(s),
	Ladder<std::complex<double> >(set_spuc(t,flux),"ladder-freeflux"),
	t_(t),
	flux_(flux)
{
	if(status_==2){
		init_fermionic();

		system_info_.text("LadderFreeFlux :");
		system_info_.item("Each color has the same Hamiltonian.");
		system_info_.item("Free hopping term");
		system_info_.item("Free flux per plaquette ");
		filename_ += "-t";
		for(unsigned int i(0);i<t_.size();i++){
			filename_ += ((t_(i)>0)?"+":"")+my::tostring(t_(i));
		}
		filename_ += "-phi";
		for(unsigned int i(0);i<flux_.size()-1;i++){
			filename_ += ((flux_(i)>0)?"+":"")+my::tostring(flux_(i));
		}
	}
}

/*{method needed for running*/
void LadderFreeFlux::compute_H(){
	H_.set(n_,n_,0);
	Matrix<int> nb;
	std::cerr<<__PRETTY_FUNCTION__<<" : undefined"<<std::endl;
	//for(unsigned int i(0);i<n_;i++){
		//nb = get_neighbourg(i);
		//if(i%spuc_){
			//if(i%2){
				//H_(i,nb(0,0)) = my::chop(std::polar(nb(0,1)*t_(k++),flux_(f++)*M_PI));
			//} else {
				//H_(i,nb(0,0)) = nb(0,1)*t_(k++);
				//H_(i,nb(1,0)) = nb(1,1)*t_(k++);
			//}
		//} else {
			///*!decoupled chains in the limit J⊥(0)=J_(1)->0, therefore, in
			// * that case the variational parameter is t⊥ (t‖ is set to 1),
			// * otherwise the inverse is done*/
			//if(J_(0)>J_(1)){
				//H_(i,nb(0,0)) = nb(0,1);
				//H_(i,nb(1,0)) = nb(1,1)*t_(k++);
			//} else {
				//H_(i,nb(0,0)) = nb(0,1)*t_(k++);
				//H_(i,nb(1,0)) = nb(1,1);
			//}
		//}
		//k = k%t_.size();
		//f = f%flux_.size();
	//}
	H_ += H_.conjugate_transpose();
}

void LadderFreeFlux::create(){
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

void LadderFreeFlux::save_param(IOFiles& w) const {
	if(w.is_binary()){
		Vector<double> param(t_.size()+flux_.size());
		std::string s("t=(");
		for(unsigned int i(0);i<t_.size()-1;i++){
			s += my::tostring(t_(i))+",";
			param(i) = t_(i);
		}
		param(t_.size()-1) = t_.back();
		s += my::tostring(t_.back())+") : " + RST::math("\\phi/\\pi")+"=(";
		for(unsigned int i(0);i<flux_.size()-1;i++){
			s += my::tostring(flux_(i))+",";
			param(i+t_.size()) = flux_(i);
		}
		param.back() = flux_.back();
		s += my::tostring(flux_.back())+")";
		w.add_to_header()->title(s,'<');
		w<<param;
		w.add_to_header()->add(system_info_.get());
	} else { w<<t_<<" "<<flux_<<" "; }
}

unsigned int LadderFreeFlux::set_spuc(Vector<double> const& t, Vector<double> const& flux){
	if((t.size()+1)%flux.size()){
		std::cerr<<__PRETTY_FUNCTION__<<" : incoherent t and flux sizes : "<<t.size()<<" "<<flux.size()<<std::endl;
		return 1;
	} else {
		switch(t.size()){
			case 2: { return 2; } break;
			case 5: { return 4; } break;
			case 8: { return 6; } break;
			case 11:{ return 8; } break;
			default:{
						std::cerr<<__PRETTY_FUNCTION__<<" : invalid t size : "<<t.size()<<std::endl;
						return 1;
					}
		}
	}
}
/*}*/

/*{method needed for checking*/
void LadderFreeFlux::display_results(){
	compute_H();
	std::string t(my::tostring(t_));
	std::string flux(my::tostring(flux_));
	draw_lattice(true,true,true,"-d:t "+t+" -d:phi "+flux, "t=("+t+")"+RST::math("\\phi")+"=("+flux+")");
}

void LadderFreeFlux::check(){
	info_ = "";
	path_ = "";
	dir_  = "./";
	filename_ ="ladder-freeflux";
	display_results();

	//plot_band_structure();
}
/*}*/

/*{method needed for analysing*/
std::string LadderFreeFlux::extract_level_6(){
	(*data_write_)<<N_<<" "<<m_<<" "<<n_<<" "<<bc_<<" "<<asin(J_(1))<<" "<<obs_[0][0]<<IOFiles::endl;

	display_results();

	save_param(*jd_write_);
	save(*jd_write_);

	return filename_;
}
/*}*/
