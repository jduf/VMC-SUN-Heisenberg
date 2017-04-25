#include "Observable.hpp"

/*constructors and destructor*/
/*{*/
Observable::Observable(std::string const& name, unsigned int const& type, unsigned int const& nval, Matrix<int> const& links, unsigned int const& modulo, unsigned int const& B, unsigned int const& b, bool const& conv):
	name_(name),
	type_(type),
	modulo_(modulo),
	nval_(nval),
	links_(links),
	val_(NULL)
{ set(B,b,conv); }

Observable::Observable(std::string const& name, unsigned int const& type, unsigned int const& nval, unsigned int const& nlinks, unsigned int const& modulo, unsigned int const& B, unsigned int const& b, bool const& conv):
	name_(name),
	type_(type),
	modulo_(modulo),
	nval_(nval),
	links_(nlinks,3),
	val_(NULL)
{ set(B,b,conv); }

Observable::Observable(IOFiles& r):
	name_(r.read<std::string>()),
	type_(r.read<unsigned int>()),
	modulo_(r.read<unsigned int>()),
	nval_(r.read<unsigned int>()),
	links_(r),
	val_(new Data<double>[nval_])
{
	for(unsigned int i(0);i<nval_;i++){ val_[i] = std::move(Data<double>(r)); }
}

Observable::Observable(Observable const& obs):
	name_(obs.name_),
	type_(obs.type_),
	modulo_(obs.modulo_),
	nval_(obs.nval_),
	links_(obs.links_),
	val_(obs.val_?new Data<double>[nval_]:NULL)
{
	for(unsigned int i(0);i<nval_;i++){ val_[i] = obs.val_[i]; }
}

Observable::Observable(Observable&& obs):
	name_(std::move(obs.name_)),
	type_(obs.type_),
	modulo_(obs.modulo_),
	nval_(obs.nval_),
	links_(std::move(obs.links_)),
	val_(obs.val_)
{
	obs.val_ = NULL;
	obs.nval_= 0;
}

void Observable::set(unsigned int const& B, unsigned int const& b, bool const& conv){
	if(!modulo_){
		if(links_.row()%nval_){
			nval_ = 0;
			std::cerr<<__PRETTY_FUNCTION__<<" : incoherent number : nval="<<nval_<<" nlinks="<<links_.row()<<" for '"<<name_<<"'"<<std::endl;
		} else { modulo_ = links_.row()/nval_; }
	}
	if(nval_){
		val_ = new Data<double>[nval_];
		for(unsigned int i(0);i<nval_;i++){ val_[i].set(B,b,conv); }
	}
}

Observable::~Observable(){
	if(val_){ delete[] val_; }
}
/*}*/

/*assignement operator*/
/*{*/
Observable& Observable::operator=(Observable obs){
	swap_to_assign(*this,obs);
	return (*this);
}

void Observable::swap_to_assign(Observable& obs1, Observable& obs2){
	std::swap(obs1.modulo_,obs2.modulo_);
	std::swap(obs1.nval_,obs2.nval_);
	std::swap(obs1.val_,obs2.val_);
	std::swap(obs1.links_,obs2.links_);
}
/*}*/

/*handles class attributes*/
/*{*/
void Observable::reset(){
	for(unsigned int i(0);i<nval_;i++){ val_[i].set(); }
}

void Observable::set_x(double const& val){
	for(unsigned int i(0);i<nval_;i++){ val_[i].set_x(val); }
}

void Observable::add(unsigned int const& i, double const& val){
	val_[i].add(val/modulo_);
}

void Observable::add_sample(){
	for(unsigned int i(0);i<nval_;i++){ val_[i].add_sample(); }
}

void Observable::merge(Observable& obs){
	if(type_ == 1234){ 
		type_ = 0; 
		std::cerr<<__PRETTY_FUNCTION__<<" : redefine the type_ of the Observable 'Energy per site'"<<std::endl;
		/*! This can is a well defined behaviour when VMCMinimization tries to
		 * measure the Obervable 'Bond energy' for a new sample. This is
		 * because for a new sample created using Minimization::s_, the
		 * Observable 'Energy per site' will have its type_ changed to 1234 and
		 * when merging with another mesure of the same sample with type 1234,
		 * this method would merge two Observables with the type :
		 * type_ = obs.type_ = 1234
		 *
		 * To keep the type_ of the Obervable 'Energy per site' equal to 0, it
		 * is redefined here.*/
	}
	/*see Observable::combine_measurement for explanation on the 1234 value*/
	if(type_ == (obs.type_)%1234 && nval_ == obs.nval_){
		for(unsigned int i(0);i<nval_;i++){ val_[i].merge(obs.val_[i]); }
	} else {
		std::cerr<<__PRETTY_FUNCTION__<<" : inconsistent type or size : "
			<<name_    <<" ("<<type_    <<","<<nval_<<") != "
			<<obs.name_<<" ("<<obs.type_<<","<<obs.nval_<<")"<<std::endl; 
	}
}

void Observable::delete_binning(){
	for(unsigned int i(0);i<nval_;i++){ val_[i].delete_binning(); }
}

void Observable::complete_analysis(double const& convergence_criterion){
	for(unsigned int i(0);i<nval_;i++){ val_[i].complete_analysis(convergence_criterion); }
}

void Observable::combine_measurement(bool const& combine){
	/*!If the measure of a type 1234 is never defined, it will never be
	 * measured. This is the behaviour needed to compute the energy and the
	 * bond energy at the same time. The energy shouldn't be computed alone but
	 * should be computed when the bond energy (with type 1) is computed*/
	if(!type_){
		if(combine){ type_ = 1234; }
		else { type_ = 0; }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : should only be used with the Observable 'Energy per site' that can be measured with the Observable 'Bond energy'"<<std::endl; }
}
/*}*/

/*reads/writes in IOFiles methods and print*/
/*{*/
IOFiles& operator<<(IOFiles& w, Observable const& obs){
	if(w.is_binary()){ obs.write(w); }
	else { w.stream()<<obs; }
	return w;
}

IOFiles& operator>>(IOFiles& r, Observable& obs){
	if(r.is_binary()){ obs = Observable(r); }
	else { obs = Observable(r); }
	return r;
}

std::ostream& operator<<(std::ostream& flux, Observable const& obs){
	flux<<obs.get_name()<<std::endl;
	for(unsigned int i(0);i<obs.nval();i++){ flux<<obs[i]<<std::endl; }
	return flux;
}

void Observable::write(IOFiles& w) const {
	w<<name_<<type_<<modulo_<<nval_<<links_;
	for(unsigned int i(0);i<nval_;i++){ w<<val_[i]; }
}

void Observable::print(bool const& all) const {
	std::cout<<name_<<std::endl;
	if(all){
		if(type_ == 0){ std::cout<<val_[0]<<std::endl; }
		else {
			for(unsigned int i(0);i<links_.row();i++){
				std::cout<<links_(i,0)<<" "<<links_(i,1)<<" "<<links_(i,2)<<" "<<val_[links_(i,2)]<<std::endl;
			}
		}
	} else { std::cout<<links_<<std::endl; }
}
/*}*/
