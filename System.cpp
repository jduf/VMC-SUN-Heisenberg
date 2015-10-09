#include "System.hpp"

/*constructors and destructor*/
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
	unsigned int n_obs_(r.read<unsigned int>());
	for(unsigned int i(0);i<n_obs_;i++){
		obs_.push_back(Observable(r));
	}
}

System::System(System const& s):
	ref_(s.ref_),
	N_(s.N_), 
	m_(s.m_),
	n_(s.n_),
	bc_(s.bc_),
	M_(s.M_),
	J_(s.J_),
	status_(s.status_),
	obs_(s.obs_)
{}


/*}*/

/*handles class attributes*/
/*{*/
void System::set(Vector<double> const& J, std::vector<Observable> const& obs, unsigned int const& which_observables){
	//E_.set(50,5,false);
	J_ = J; 
	for(unsigned int i(0);i<which_observables;i++){
		obs_.push_back(obs[i]);
	}
}

void System::clear_measurments(){ 
	//E_.set(); 
	for(unsigned int i(0);i<obs_.size();i++){ obs_[i].set(); }
}

bool System::check_conv(double const& convergence_criterion){
	(void)(convergence_criterion);
	//E_.complete_analysis(convergence_criterion);
	//return E_.get_conv();
	return 0;
}

void System::complete_analysis(double const& convergence_criterion){ 
	//E_.complete_analysis(convergence_criterion); 
	for(unsigned int i(0);i<obs_.size();i++){
		obs_[i].complete_analysis(convergence_criterion); 
	}
}

void System::merge(System* s){
	//E_.merge(s->E_);
	for(unsigned int i(0);i<obs_.size();i++){
		obs_[i].merge(s->obs_[i]); 
	}
}

void System::delete_binning(){
	//E_.delete_binning();
	for(unsigned int i(0);i<obs_.size();i++){
		obs_[i].delete_binning();
	}
}
/*}*/

/*write in IOFiles methods*/
/*{*/
void System::write(IOFiles& w) const {
	if(w.is_binary()){
		unsigned int nobs(obs_.size());
		w<<ref_<<N_<<m_<<n_<<bc_<<M_<<J_<<status_<<nobs;
		for(unsigned int i(0);i<nobs;i++){ w<<obs_[i]; }
	} else {
		w<<N_<<" "<<m_<<" "<<n_<<" "<<bc_<<" "<<M_<<IOFiles::endl;
	}
}

void System::save_input(IOFiles& w) const {
	RST rst;
	rst.title("Input",'+');
	w.add_header()->add(rst.get());

	w.write("ref (type of wavefunction)",ref_);
	w.write("N (N of SU(N))",N_);
	w.write("m (# of particles per site)",m_);
	w.write("n (# of site)",n_);
	w.write("bc (boundary condition)",bc_);
	w.write("M (# of particles of each color, "+my::tostring(M_(0))+")",M_);
	w.write("J (energy of each bond)",J_);
}

void System::save_output(IOFiles& w) const {
	RST rst;
	rst.title("Results",'+');
	w.add_header()->add(rst.get());

	w.write("status",status_);
	unsigned int nobs(obs_.size());
	w.write("number of types of correlations (saved as well)",nobs);
	for(unsigned int i(0);i<corr_types_.size();i++){ w<<obs_[i]; }
}
/*}*/

Vector<unsigned int> System::set_ref(Parseur& P){
	std::string const& wf(P.get<std::string>("wf"));
	Vector<unsigned int> ref(3,6);
	if( wf == "chainfermi" ){
		ref(0) = 1;
		ref(1) = 1;
		ref(2) = 0;

		std::vector<double> Jp(1,1);
		P.set("Jp",Jp);
	}
	if( wf == "chainpolymerized" ){
		ref(0) = 1;
		ref(1) = 1;
		ref(2) = 1;

		std::vector<double> Jp(1,1);
		P.set("Jp",Jp);
	}

	if( wf == "ladderfermi"){
		ref(0) = 2;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "ladderfree"){
		ref(0) = 2;
		ref(1) = 1;
		ref(2) = 4;
		std::vector<double> Jp(2);
		double theta(P.get<double>("theta"));
		Jp[0] = cos(theta); //legs  (J‖)
		Jp[1] = sin(theta); //rungs (J⊥)
		P.set("Jp",Jp);
	}

	if( wf == "trianglefermi" ){
		ref(0) = 3;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "trianglemu" ){
		ref(0) = 3;
		ref(1) = 1;
		ref(2) = 1;
	}
	if( wf == "trianglephi" ){
		ref(0) = 3;
		ref(1) = 2;
		ref(2) = 2;
	}
	if( wf == "trianglejastrow" ){
		ref(0) = 3;
		ref(1) = 0;
	}

	if( wf == "squarefermi" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "squaremu" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 1;
	}
	//if( wf == "squarefreereal" ){
	//ref(0) = 4;
	//ref(1) = 1;
	//ref(2) = 3;
	//C_.set("t",param?param->range(0,3):C->get<std::vector<double> >("t"));
	//C_.set("mu",param?param->range(3,5):C->get<std::vector<double> >("mu"));
	//}
	if( wf == "squareacsl" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 3;
		std::vector<double> Jp(1,1);
		P.set("Jp",Jp);
	}
	if( wf == "squarefreecomplex" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 4;
		std::vector<double> Jp(1,1);
		P.set("Jp",Jp);
	}
	if( wf == "squarepiflux" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 2;
		std::vector<double> Jp(1,1);
		P.set("Jp",Jp);
	}
	if( wf == "squarephi" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 3;
	}
	if( wf == "squarejastrow" ){
		ref(0) = 4;
		ref(1) = 0;
	}

	if( wf == "kagomefermi" ){
		ref(0) = 5;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "kagomedirac" ){
		ref(0) = 5;
		ref(1) = 1;
		ref(2) = 1;
	}
	if( wf == "kagomevbc" ){
		ref(0) = 5;
		ref(1) = 2;
		ref(2) = 2;
	}

	if( wf == "honeycomb0pp" ){
		ref(0) = 6;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "honeycombsu4" ){
		ref(0) = 6;
		ref(1) = 1;
		ref(2) = 0;
	}
	return ref;
}
