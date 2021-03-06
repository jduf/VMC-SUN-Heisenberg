#include "System.hpp"

/*constructors*/
/*{*/
System::System(Parseur& P):
	ref_(complete_system_info(P)),
	N_(P.get<unsigned int>("N")),
	m_(P.get<unsigned int>("m")),
	n_(P.get<unsigned int>("n")),
	bc_(P.get<int>("bc")),
	M_(P.get<std::vector<unsigned int> >("M")),
	J_(P.get<std::vector<double> >("J")),
	status_(5)
{
	if(M_.sum() != m_*n_ || m_>N_){ std::cerr<<__PRETTY_FUNCTION__<<" : N, M, m and n are incompatible"<<std::endl; }
	else if(bc_ != -1 && bc_ != 0 && bc_ != 1){ std::cerr<<__PRETTY_FUNCTION__<<" : unknown boundary condition"<<std::endl; }
	else { status_ = 4; }
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
	unsigned int nobs(r.read<unsigned int>());
	for(unsigned int i(0);i<nobs;i++){ obs_.push_back(Observable(r)); }
}
/*}*/

/*handles class attributes*/
/*{*/
bool System::try_other_geometry(Vector<unsigned int> const& ref) const {
	if(ref_(0) == ref(0) && ref_(1) == ref(1) && ref_(2) == ref(2)){
		if(ref_(0) == 4){
			ref_(3) = ref(3);
			std::cerr<<__PRETTY_FUNCTION__<<std::endl;
			return ref_(3)<3;
		}
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : the wavefunction should not have changed"<<std::endl; }
	std::cerr<<__PRETTY_FUNCTION__<<" : no solution found"<<std::endl;
	return false;
}

bool System::check_conv(){
	/*check only the energy*/
	obs_[0][0].complete_analysis(1e-5);
	return obs_[0][0].get_conv();
}

void System::complete_analysis(){
	for(auto& o:obs_){ o.complete_analysis(); }
}

void System::merge(System* const s){
	/*this commented bit should be used instead of the last bit when one is
	 * trying to merge two independant simulations where the order of the
	 * observable might be different*/
	//unsigned int match[Observable::number_of_observables_defined] = { 0 };
	//for(unsigned int i(0);i<obs_.size();i++){ match[obs_[i].get_type()%1234]=i+1; }
	//for(unsigned int i(0);i<s->obs_.size();i++){
	//unsigned int o(match[s->obs_[i].get_type()%1234]);
	//if(o){ obs_[o-1].merge(s->obs_[i]); }
	//else { obs_.push_back(s->obs_[i]); }
	//}

	unsigned int i(0);
	while(i<obs_.size())   { obs_[i].merge(s->obs_[i]);  i++; }
	while(i<s->obs_.size()){ obs_.push_back(s->obs_[i]); i++; }
}

void System::delete_binning(){
	for(auto& o:obs_){ o.delete_binning(); }
}
/*}*/

/*write in IOFiles methods and print*/
/*{*/
void System::write(IOFiles& w) const {
	if(w.is_binary()){
		w<<ref_<<N_<<m_<<n_<<bc_<<M_<<J_<<status_<<(unsigned int)(obs_.size());
		for(auto& o:obs_){ w<<o; }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : not adapted to binary file :"<<w.get_filename()<<std::endl; }
}

void System::save(IOFiles& w) const {
	if(w.is_binary()){
		save_input(w);
		save_output(w);
	} else { w<<N_<<" "<<m_<<" "<<n_<<" "<<bc_<<" "<<obs_[0][0]<<" "<<ref_<<IOFiles::endl; }
}

void System::save_input(IOFiles& w) const {
	RST rst;
	rst.title("Input",'+');
	w.add_to_header()->add(rst.get());

	w.write("ref (geometry "+my::tostring(ref_(3))+")",ref_);
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
	w.add_to_header()->add(rst.get());

	w.write("status",status_);
	unsigned int nobs(obs_.size());
	if(nobs){
		w.write("#obs (E="+my::tostring(obs_[0][0].get_x())+", dE="+my::tostring(obs_[0][0].get_dx())+")",nobs);
		for(auto& o:obs_){ w<<o; }
	} else { w.write("#obs",nobs); }
}

void System::print(bool const& all) const {
	if(all){
		std::cout<<"SU("<<N_<<") m="<<m_<<" n="<<n_<<" BC="<<bc_<<" nobs="<<obs_.size()<<" ref="<<ref_<<std::endl;
		for(auto& o:obs_){ std::cout<<std::endl<<o; }
	} else { std::cout<<obs_[0][0]<<std::endl; }
	std::cout<<RST::dash_line_<<std::endl;
}
/*}*/

Vector<unsigned int> System::complete_system_info(Parseur& P){
	std::string const& wf(P.get<std::string>("wf"));
	Vector<unsigned int> ref(5,0);
	ref(4) = 2;//!to force the creation of the cluster
	if( wf == "chain-fermi" ){
		ref(0) = 1;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "chain-free" ){
		ref(0) = 1;
		ref(1) = 1;
		ref(2) = 1;
	}
	if( wf == "chain-polymerized" ){
		ref(0) = 1;
		ref(1) = 1;
		ref(2) = 2;
	}
	if( wf == "chain-sas" ){
		ref(0) = 1;
		ref(1) = 1;
		ref(2) = 3;
		std::vector<double> J(2);
		J[0] = P.get<double>("J0");
		J[1] = P.get<double>("J1");
		P.set("J",J);
	}

	if( wf == "ladder-fermi" ){
		ref(0) = 2;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "ladder-free" ){
		ref(0) = 2;
		ref(1) = 1;
		ref(2) = 1;
	}
	if( wf == "ladder-dimerA" ){
		ref(0) = 2;
		ref(1) = 1;
		ref(2) = 2;
	}
	if( wf == "ladder-dimerB" ){
		ref(0) = 2;
		ref(1) = 1;
		ref(2) = 3;
	}
	if( wf == "ladder-squareplaquetteA" ){
		ref(0) = 2;
		ref(1) = 1;
		ref(2) = 4;
	}
	if( wf == "ladder-squareplaquetteB" ){
		ref(0) = 2;
		ref(1) = 1;
		ref(2) = 5;
	}
	if( wf == "ladder-squareplaquetteC" ){
		ref(0) = 2;
		ref(1) = 1;
		ref(2) = 6;
	}
	if( wf == "ladder-rectangularplaquetteA" ){
		ref(0) = 2;
		ref(1) = 1;
		ref(2) = 7;
	}
	if( wf == "ladder-rectangularplaquetteB" ){
		ref(0) = 2;
		ref(1) = 1;
		ref(2) = 8;
	}
	if( wf == "ladder-rectangularplaquetteC" ){
		ref(0) = 2;
		ref(1) = 1;
		ref(2) = 9;
	}
	if( wf == "ladder-rectangularplaquetteD" ){
		ref(0) = 2;
		ref(1) = 1;
		ref(2) = 10;
	}

	if( wf == "triangle-fermi" ){
		ref(0) = 3;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "triangle-free" ){
		ref(0) = 3;
		ref(1) = 1;
		ref(2) = 1;
	}
	if( wf == "triangle-plaquette" ){
		ref(0) = 3;
		ref(1) = 1;
		ref(2) = 2;
	}
	if( wf == "triangle-mu" ){
		ref(0) = 3;
		ref(1) = 1;
		ref(2) = 3;
	}
	if( wf == "triangle-T3x2" ){
		ref(0) = 3;
		ref(1) = 1;
		ref(2) = 4;
	}
	if( wf == "triangle-alternatingplaquette" ){
		ref(0) = 3;
		ref(1) = 1;
		ref(2) = 5;
	}
	if( wf == "triangle-chiral" ){
		ref(0) = 3;
		ref(1) = 2;
		ref(2) = 1;
	}
	if( wf == "triangle-phi" ){
		ref(0) = 3;
		ref(1) = 2;
		ref(2) = 2;
	}
	if( wf == "triangle-chiralSG" ){
		ref(0) = 3;
		ref(1) = 2;
		ref(2) = 3;
	}

	if( wf == "square-fermi" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "square-free" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 1;
	}
	if( wf == "square-mu" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 2;
	}
	if( wf == "square-dimerizedbar" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 3;
	}
	if( wf == "square-T2x2" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 4;
	}
	if( wf == "square-T3x2" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 5;
	}
	if( wf == "square-T3x3" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 6;
	}
	if( wf == "square-T4x2" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 7;
	}
	if( wf == "square-T4x3" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 8;
	}
	if( wf == "square-T4x4" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 9;
	}
	if( wf == "square-ladder" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 10;
	}
	if( wf == "square-VCS" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 11;
	}
	if( wf == "square-mu-T2x1" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 12;
	}
	if( wf == "square-mu-T2x2" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 13;
	}
	if( wf == "square-facingdimers" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 14;
	}
	if( wf == "square-alternatingdimers" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 15;
	}
	if( wf == "square-k2mupi" ){
		ref(0) = 4;
		ref(1) = 1;
		ref(2) = 16;
	}
	if( wf == "square-freeflux" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 1;
	}
	if( wf == "square-piflux" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 2;
	}
	if( wf == "square-chiral" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 3;
	}
	if( wf == "square-box6" ){
		ref(0) = 4;
		ref(1) = 2;
		ref(2) = 4;
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
	if( wf == "kagome-free" ){
		ref(0) = 5;
		ref(1) = 1;
		ref(2) = 1;
	}
	if( wf == "kagome-plaquette3A" ){
		ref(0) = 5;
		ref(1) = 1;
		ref(2) = 2;
	}
	if( wf == "kagome-plaquette3B" ){
		ref(0) = 5;
		ref(1) = 1;
		ref(2) = 3;
	}
	if( wf == "kagome-plaquette6A" ){
		ref(0) = 5;
		ref(1) = 1;
		ref(2) = 4;
	}
	if( wf == "kagome-plaquette6B" ){
		ref(0) = 5;
		ref(1) = 1;
		ref(2) = 5;
	}
	if( wf == "kagome-chiral" ){
		ref(0) = 5;
		ref(1) = 2;
		ref(2) = 1;
	}
	if( wf == "kagome-vbc" ){
		ref(0) = 5;
		ref(1) = 2;
		ref(2) = 2;
	}
	if( wf == "kagome-pihalftriangle" ){
		ref(0) = 5;
		ref(1) = 2;
		ref(2) = 3;
	}
	if( wf == "kagome-chiralB" ){
		ref(0) = 5;
		ref(1) = 2;
		ref(2) = 4;
	}

	if( wf == "honeycomb-fermi" ){
		ref(0) = 6;
		ref(1) = 1;
		ref(2) = 0;
	}
	if( wf == "honeycomb-free" ){
		ref(0) = 6;
		ref(1) = 1;
		ref(2) = 1;
	}
	if( wf == "honeycomb-0pp" ){
		ref(0) = 6;
		ref(1) = 1;
		ref(2) = 3;
	}
	if( wf == "honeycomb-p00" ){
		ref(0) = 6;
		ref(1) = 1;
		ref(2) = 4;
	}
	if( wf == "honeycomb-piflux" ){
		ref(0) = 6;
		ref(1) = 1;
		ref(2) = 5;
	}
	if( wf == "honeycomb-chiral" ){
		ref(0) = 6;
		ref(1) = 2;
		ref(2) = 1;
	}

	switch(ref(0)){
		case 2:
			{
				std::vector<double> J(2);
				double theta(P.get<double>("theta"));
				J[0] = cos(theta); //legs  (J‖)
				J[1] = sin(theta); //rungs (J⊥)
				P.set("J",J);
			}break;
		case 4:
			{
				unsigned int n(P.get<unsigned int>("n"));
				unsigned int i;
				if(P.find("geometry",i)){ ref(3) = P.get<unsigned int>(i); }
				else {
					ref(3)=1;
					if(my::are_equal(sqrt(n),floor(sqrt(n)))){
						for(unsigned int p(0);p<=sqrt(n);p++){
							for(unsigned int q(1);q<p+1;q++){
								if(p*p+q*q==n){
									std::cerr<<__PRETTY_FUNCTION__<<" : tilted geometry chosen by default (to change -u:geometry 0)"<<std::endl;
									q = p = n;
								}
							}
						}
					}
				}
				if(P.find("Jp")){
					std::vector<double> J(2);
					double Jp(P.get<double>("Jp"));
					J[0] = 1; //horizontal links
					J[1] = Jp;//vertical links
					P.set("J",J);
				}
			}break;
		case 5:
			{
				unsigned int n(P.get<unsigned int>("n"));
				if(my::are_equal(sqrt(n/9.0),floor(sqrt(n/9.0)))){ ref(3) = 1; }
			}break;
		case 6:
			{
				unsigned int n(P.get<unsigned int>("n"));
				ref(3) = 2;
				if(my::are_equal(sqrt(n/2.0),floor(sqrt(n/2.0)))){ ref(3) = 0; }
				if(my::are_equal(sqrt(n/6.0),floor(sqrt(n/6.0)))){ ref(3) = 1; }
			}break;
	}
	if(!P.find("J")){ P.set("J",std::vector<double>(1,1)); }
	if(!P.find("M")){ P.set("M",std::vector<unsigned int>(P.get<unsigned int>("N"),P.get<unsigned int>("n")*P.get<unsigned int>("m")/P.get<unsigned int>("N"))); }
	return ref;
}
