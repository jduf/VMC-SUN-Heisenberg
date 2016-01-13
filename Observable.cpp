#include "Observable.hpp"

/*constructors and destructor*/
/*{*/
Observable::Observable(std::string const& name, unsigned int const& type, unsigned int const& nval, Matrix<int> const& links, unsigned int const& B, unsigned int const& b, bool const& conv):
	name_(name),
	type_(type),
	modulo_(0),
	nval_(nval),
	links_(links),
	val_(NULL)
{ set(B,b,conv); }

Observable::Observable(std::string const& name, unsigned int const& type, unsigned int const& nval, unsigned int const& nlinks, unsigned int const& B, unsigned int const& b, bool const& conv):
	name_(name),
	type_(type),
	modulo_(0),
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
	if(links_.row()%nval_){ std::cerr<<__PRETTY_FUNCTION__<<" : incoherent number"<<std::endl; }
	else {
		modulo_ = links_.row()/nval_;
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
void Observable::set_x(double const& val){
	for(unsigned int i(0);i<nval_;i++){
		val_[i].set_x(val);
	}
}

void Observable::add(unsigned int const& i, double const& val){
	val_[i].add(val/modulo_);
}

void Observable::add_sample(){
	for(unsigned int i(0);i<nval_;i++){ val_[i].add_sample(); }
}

void Observable::merge(Observable& obs){
	if(nval_ == obs.nval_){
		for(unsigned int i(0);i<nval_;i++){ val_[i].merge(obs.val_[i]); }
	} else { std::cerr<<__PRETTY_FUNCTION__<<" : inconsistent size"<<std::endl; }
}

void Observable::delete_binning(){
	for(unsigned int i(0);i<nval_;i++){ val_[i].delete_binning(); }
}

void Observable::complete_analysis(double const& convergence_criterion){
	for(unsigned int i(0);i<nval_;i++){
		val_[i].complete_analysis(convergence_criterion);
	}
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

void Observable::print() const {
	std::cout<<name_<<std::endl;
	if(nval_>1){
		for(unsigned int i(0);i<links_.row();i++){
			std::cout<<links_(i,0)<<" "<<links_(i,1)<<" "<<links_(i,2)<<" "<<val_[links_(i,2)]<<std::endl;
		}
	} else { std::cout<<links_<<std::endl; }
}
/*}*/
