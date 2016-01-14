#include "System.hpp"

/*constructors*/
/*{*/
System::System(Parseur& P):
	ref_(set_ref(P)),
	N_(P.get<unsigned int>("N")),
	m_(P.get<unsigned int>("m")),
	n_(P.get<unsigned int>("n")),
	bc_(P.get<int>("bc")),
	M_(P.get<std::vector<unsigned int> >("M")),
	J_(P.get<std::vector<double> >("Jp")),
	status_(5)
{
	if(M_.sum() != m_*n_ || m_>N_){ std::cerr<<__PRETTY_FUNCTION__<<" : Bad initialization"<<std::endl; }
	else{ status_--; }
	if(bc_ != -1 && bc_ != 0 && bc_ != 1){ std::cerr<<__PRETTY_FUNCTION__<<" : unknown boundary condition"<<std::endl; }
	else { status_--; }
}

System::System(IOFiles& r):
	ref_(r),
	N_(r.read<unsigned int>()),
	m_(r.read<unsigned int>()),
	n_(r.read<unsigned int>()),
	bc_(r.read<int>()),
	M_(r),
	J_(r),
	status_(r.read<unsigned int>())
{
	int nobs(r.read<int>());
	for(int i(0);i<nobs;i++){ obs_.push_back(Observable(r)); }
}
/*}*/

/*handles class attributes*/
/*{*/
void System::set_obs(std::vector<Observable> const& obs, int const& nobs){
	if(nobs<0){ obs_ = obs; }
	else {
		obs_.clear();
		for(int i(0);i<nobs;i++){ obs_.push_back(obs[i]); }
	}
}

void System::clear_obs(unsigned int const& from){
	for(unsigned int i(from);i<obs_.size();i++){ obs_.pop_back(); }
}

void System::reset_obs(){
	for(unsigned int i(0);i<obs_.size();i++){ obs_[i].reset(); }
}

bool System::check_conv(double const& convergence_criterion){
	/*check only the energy*/
	obs_[0][0].complete_analysis(convergence_criterion);
	return obs_[0][0].get_conv();
}

void System::complete_analysis(double const& convergence_criterion){
	for(unsigned int i(0);i<obs_.size();i++){ obs_[i].complete_analysis(convergence_criterion); }
}

void System::merge(System* const s){
	unsigned int i(0);
	while(i<obs_.size()){ obs_[i].merge(s->obs_[i]); i++; }
	while(i<s->obs_.size()){ obs_.push_back(s->obs_[i]); i++; }
}

void System::delete_binning(){
	for(unsigned int i(0);i<obs_.size();i++){ obs_[i].delete_binning(); }
}
/*}*/

/*write in IOFiles methods and print*/
/*{*/
void System::write(IOFiles& w) const {
	int nobs(obs_.size());
	w<<ref_<<N_<<m_<<n_<<bc_<<M_<<J_<<status_<<nobs;
	for(int i(0);i<nobs;i++){ w<<obs_[i]; }
}

void System::save(IOFiles& w) const {
	if(w.is_binary()){
		save_input(w);
		save_output(w);
	} else { w<<N_<<" "<<m_<<" "<<n_<<" "<<bc_<<" "<<M_<<" "<<obs_[0][0]<<IOFiles::endl; }
}

void System::save_input(IOFiles& w) const {
	RST rst;
	rst.title("Input",'+');
	w.add_header()->add(rst.get());

	w.write("ref (wavefunction)",ref_);
	w.write("SU",N_);
	w.write("m (#particles per site)",m_);
	w.write("n (#site)",n_);
	w.write("bc (boundary condition)",bc_);
	w.write("M (#particles per color, "+my::tostring(M_(0))+")",M_);
	w.write("J (bond strength)",J_);
}

void System::save_output(IOFiles& w) const {
	RST rst;
	rst.title("Results",'+');
	w.add_header()->add(rst.get());

	w.write("status",status_);
	int nobs(obs_.size());
	w.write("#obs (E="+my::tostring(obs_[0][0].get_x())+", dE="+my::tostring(obs_[0][0].get_dx())+")",nobs);
	for(int i(0);i<nobs;i++){ w<<obs_[i]; }
}

void System::print(unsigned int nobs) const {
	std::cout<<std::string(35,'-')<<std::endl;
	std::cout<<"SU("<<N_<<") m="<<m_<<" n="<<n_<<" BC="<<bc_<<" nobs="<<obs_.size()<<std::endl<<std::endl;
	if(nobs>obs_.size()){ nobs = obs_.size(); }
	for(unsigned int i(0);i<nobs;i++){ std::cout<<obs_[i]<<std::endl; }
}
/*}*/

Vector<unsigned int> System::set_ref(Parseur& P){
	std::string const& wf(P.get<std::string>("wf"));
	Vector<unsigned int> ref(3,6);
	if( wf == "chain-fermi" ){
		ref(0) = 1;
		ref(1) = 1;
		ref(2) = 0;

		std::vector<double> Jp(1,1);
		P.set("Jp",Jp);
	}
	if( wf == "chain-polymerized" ){
		ref(0) = 1;
		ref(1) = 1;
		ref(2) = 1;

		std::vector<double> Jp(1,1);
		P.set("Jp",Jp);
	}
	if( wf == "chain-free" ){
		ref(0) = 1;
		ref(1) = 1;
		ref(2) = 2;
	}

	if( wf == "ladder-fermi" ){
		ref(0) = 2;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "ladder-free" ){
		ref(0) = 2;
		ref(1) = 1;
		ref(2) = 4;
		std::vector<double> Jp(2);
		double theta(P.get<double>("theta"));
		Jp[0] = cos(theta); //legs  (J‖)
		Jp[1] = sin(theta); //rungs (J⊥)
		P.set("Jp",Jp);
	}
	if( wf == "ladder-freeflux" ){
		ref(0) = 2;
		ref(1) = 2;
		ref(2) = 1;
		std::vector<double> Jp(2);
		double theta(P.get<double>("theta"));
		Jp[0] = cos(theta); //legs  (J‖)
		Jp[1] = sin(theta); //rungs (J⊥)
		P.set("Jp",Jp);
	}

	if( wf == "triangle-fermi" ){
		ref(0) = 3;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "triangle-mu" ){
		ref(0) = 3;
		ref(1) = 1;
		ref(2) = 1;
	}
	if( wf == "triangle-phi" ){
		ref(0) = 3;
		ref(1) = 2;
		ref(2) = 2;
	}
	if( wf == "triangle-jastrow" ){
		ref(0) = 3;
		ref(1) = 0;
	}

	if( wf == "square-fermi" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "square-mu" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 1;
	}
	//if( wf == "square-freereal" ){
	//ref(0) = 4;
	//ref(1) = 1;
	//ref(2) = 3;
	//C_.set("t",param?param->range(0,3):C->get<std::vector<double> >("t"));
	//C_.set("mu",param?param->range(3,5):C->get<std::vector<double> >("mu"));
	//}
	if( wf == "square-acsl" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 3;
		std::vector<double> Jp(1,1);
		P.set("Jp",Jp);
	}
	if( wf == "square-freeflux" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 4;
		std::vector<double> Jp(1,1);
		P.set("Jp",Jp);
	}
	if( wf == "square-piflux" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 2;
		std::vector<double> Jp(1,1);
		P.set("Jp",Jp);
	}
	if( wf == "square-phi" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 3;
	}
	if( wf == "square-jastrow" ){
		ref(0) = 4;
		ref(1) = 0;
	}

	if( wf == "kagome-fermi" ){
		ref(0) = 5;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "kagome-dirac" ){
		ref(0) = 5;
		ref(1) = 1;
		ref(2) = 1;
	}
	if( wf == "kagome-vbc" ){
		ref(0) = 5;
		ref(1) = 2;
		ref(2) = 2;
	}

	if( wf == "honeycomb-0pp" ){
		ref(0) = 6;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "honeycomb-piflux" ){
		ref(0) = 6;
		ref(1) = 1;
		ref(2) = 1;
	}
	return ref;
}
