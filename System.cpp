#include "System.hpp"

/*constructors and destructor*/
/*{*/
System::System(Parseur& P):
	ref_(set_ref(P)),
	N_(P.get<unsigned int>("N")),
	m_(P.get<unsigned int>("m")),
	n_(P.get<unsigned int>("n")),
	bc_(P.get<int>("bc")),
	M_(P.get<std::vector<unsigned int>>("M")),
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
	status_(r.read<unsigned int>()),
	links_(r),
	E_(r),
	corr_(r),
	lr_corr_(r)
{}

System::System(System const& s, unsigned int const& status):
	ref_(s.ref_),
	N_(s.N_),
	m_(s.m_),
	n_(s.n_),
	bc_(s.bc_),
	M_(s.M_),
	J_(s.J_),
	status_(status),
	links_(s.links_),
	E_(s.E_),
	corr_(s.corr_),
	lr_corr_(s.lr_corr_)
{}
/*}*/

void System::set_bonds(System const* const s){
	J_ = s->J_; 
	links_ = s->links_; 
}

void System::set_observable(unsigned int const& which){
	E_.set(50,5,false);
	if(which>0){ corr_.set(links_.row(),50,5,false); }
	if(which>1){ lr_corr_.set(n_,50,5,false); }
}

void System::set_binning(){ 
	E_.set(); 
	corr_.set(); 
	lr_corr_.set(); 
}

void System::delete_binning(){ 
	E_.delete_binning();
	corr_.delete_binning();
	lr_corr_.delete_binning();
}

void System::write(IOFiles& w) const {
	w<<ref_<<N_<<m_<<n_<<bc_<<M_<<J_<<status_<<links_<<E_<<corr_<<lr_corr_;
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
	w.write("links",links_);
	w.write("energy per site",E_);
	w.write("correlation on links",corr_);
	w.write("long range correlation",lr_corr_);
}

Vector<unsigned int> System::set_ref(Parseur& P){
	std::string const& wf(P.get<std::string>("wf"));
	Vector<unsigned int> ref(3,6);
	if( wf == "chainfermi" ){
		ref(0) = 1;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "chainpolymerized" ){
		ref(0) = 1;
		ref(1) = 1;
		ref(2) = 1;
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
	}
	if( wf == "squarefreecomplex" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 4;
	}
	if( wf == "squarepiflux" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 2;
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

std::vector<std::string> System::names() const {
	std::vector<std::string> names;
	names.push_back("N" + my::tostring(N_));
	names.push_back("m" + my::tostring(m_));
	names.push_back("n" + my::tostring(n_));
	std::string tmp("M");
	for(unsigned int i(0);i<this->M_.size();i++){
		tmp  += "_" + my::tostring(M_(i));
	}
	names.push_back(tmp);
	switch(bc_){
		case -1:{ names.push_back("A"); }break;
		case 0: { names.push_back("O"); }break;
		case 1: { names.push_back("P"); }break;
	}
	if(my::are_equal(ref_(0),2) && my::are_equal(ref_(1),1) && my::are_equal(ref_(2),4)){
		if(J_.size()==links_.row() || J_.size()==2){
			tmp = "theta"+my::tostring(acos(J_(0)));
		} else {
			std::cerr<<__PRETTY_FUNCTION__<<" : J_ has an incoherent size"<<std::endl;
		}
	} else { 
		tmp = "J"; 
		for(unsigned int i(0);i<this->J_.size();i++){
			tmp  += (this->J_(i)>=0?"+":"") + my::tostring(this->J_(i));
		}
	}
	names.push_back(tmp);
	names.push_back(my::tostring(this->ref_(0))+my::tostring(this->ref_(1))+my::tostring(this->ref_(2)));
	return names;
}
