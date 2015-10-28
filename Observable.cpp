#include "Observable.hpp"

Observable::Observable(unsigned int const& nlinks):
	links_(nlinks,3),
	modulo_(0)
{}

Observable::Observable(unsigned int const& nlinks, unsigned int const& nval, unsigned int const& B, unsigned int const& b, bool const& conv):
	links_(nlinks,3),
	modulo_(0)
{
	set(nval,B,b,conv);
}

Observable::Observable(IOFiles& r):
	links_(r),
	val_(r)
{
	modulo_ = links_.row()/val_.size();
}

void Observable::set(unsigned int const& nval, unsigned int const& B, unsigned int const& b, bool const& conv){
	if(links_.row()%nval){ std::cerr<<__PRETTY_FUNCTION__<<" : incoherent number"<<std::endl; }
	else {
		modulo_ = links_.row()/nval;
		val_.set(nval,B,b,conv);
	}
}

IOFiles& operator<<(IOFiles& w, Observable const& obs){
	if(w.is_binary()){ obs.write(w); }
	return w;
}

IOFiles& operator>>(IOFiles& r, Observable& obs){
	if(r.is_binary()){ obs = Observable(r); }
	return r;
}

void Observable::swap_to_assign(Observable& obs1, Observable& obs2){
	std::swap(obs1.links_,obs2.links_);
	std::swap(obs1.val_,obs2.val_);
}

Observable& Observable::operator=(Observable obs){
	swap_to_assign(*this,obs);
	return (*this);
}

void Observable::write(IOFiles& w) const {
	w<<links_<<val_;
}

std::ostream& operator<<(std::ostream& flux, Observable const& obs){
	for(unsigned int i(0);i<obs.nval();i++){
		flux<<obs(i,0)<<" "<<obs(i,1)<<" "<<obs(i,2)<<" "<<obs[i]<<std::endl;
	}
	return flux;
}

void Observable::set_x(double const& val){
	for(unsigned int i(0);i<val_.size();i++){ val_[i].set_x(val); }
}

void Observable::add(unsigned int const& i, double const& val){
	val_[links_(i,2)].add(val/modulo_);
}

void Observable::add_sample(){ val_.add_sample(); }

void Observable::print() const {
	for(unsigned int i(0);i<links_.row();i++){
		std::cout<<links_(i,0)<<" "<<links_(i,1)<<" "<<links_(i,2)<<" "<<val_[links_(i,2)]<<std::endl;
	}
}

