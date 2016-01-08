#include "Observable.hpp"

/*constructors and destructor*/
/*{*/
Observable::Observable(Matrix<int> const& links, unsigned int const& nval, unsigned int const& B, unsigned int const& b, bool const& conv):
	nval_(nval),
	val_(NULL),
	links_(links)
{ set(B,b,conv); }

Observable::Observable(unsigned int const& nlinks, unsigned int const& nval, unsigned int const& B, unsigned int const& b, bool const& conv):
	nval_(nval),
	val_(NULL),
	links_(nlinks,3)
{ set(B,b,conv); }

Observable::Observable(IOFiles& r):
	nval_(r.read<unsigned int>()),
	val_(new Data<double>[nval_]),
	links_(r)
{
	if(nval_){
		modulo_ = links_.row()/nval_;
		for(unsigned int i(0);i<nval_;i++){
			val_[i] = std::move(Data<double>(r));
		}
	}
}

Observable::Observable(Observable const& obs):
	modulo_(obs.modulo_),
	nval_(obs.nval_),
	val_(obs.val_?new Data<double>[nval_]:NULL),
	links_(obs.links_)
{
	for(unsigned int i(0);i<nval_;i++){ val_[i] = obs.val_[i]; }
}

Observable::Observable(Observable&& obs):
	modulo_(obs.modulo_),
	nval_(obs.nval_),
	val_(obs.val_),
	links_(std::move(obs.links_))
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

/*read/write in IOFiles methods and print*/
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
	for(unsigned int i(0);i<obs.nval();i++){
		flux<<obs(i,0)<<" "<<obs(i,1)<<" "<<obs(i,2)<<" "<<obs[i]<<std::endl;
	}
	return flux;
}

void Observable::write(IOFiles& w) const {
	std::cout<<"check alk"<<__PRETTY_FUNCTION__<<std::endl;
	w<<links_<<val_;
}

void Observable::print() const {
	if(nval_>1){
		for(unsigned int i(0);i<links_.row();i++){
			std::cout<<links_(i,0)<<" "<<links_(i,1)<<" "<<links_(i,2)<<" "<<val_[links_(i,2)]<<std::endl;
		}
	} else { std::cout<<links_<<std::endl; }
}
/*}*/

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
