#include "Observable.hpp"

/*constructors and destructor*/
/*{*/
Observable::Observable(unsigned int const& nlinks, unsigned int const& nval, unsigned int const& B, unsigned int const& b, bool const& conv):
	links_(nlinks,3),
	modulo_(0)
{ set(nval,B,b,conv); }

Observable::Observable(IOFiles& r):
	links_(r),
	val_(r)
{
	if(val_.size()){ modulo_ = links_.row()/val_.size(); }
	else { modulo_ = 0; }
}

Observable::Observable(Matrix<int> const& links):
	links_(links),
	modulo_(0)
{}

void Observable::set(unsigned int nval, unsigned int const& B, unsigned int const& b, bool const& conv){
	if(links_.row()%nval){ std::cerr<<__PRETTY_FUNCTION__<<" : incoherent number"<<std::endl; }
	else {
		modulo_ = links_.row()/nval;
		val_.set(nval,B,b,conv);
	}
}
/*}*/

/*assignement operator*/
/*{*/
Observable& Observable::operator=(Observable obs){
	swap_to_assign(*this,obs);
	return (*this);
}

void Observable::swap_to_assign(Observable& obs1, Observable& obs2){
	std::swap(obs1.links_,obs2.links_);
	std::swap(obs1.val_,obs2.val_);
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
	else { obs = std::move(Observable(r)); }
	return r;
}

std::ostream& operator<<(std::ostream& flux, Observable const& obs){
	for(unsigned int i(0);i<obs.nval();i++){
		flux<<obs(i,0)<<" "<<obs(i,1)<<" "<<obs(i,2)<<" "<<obs[i]<<std::endl;
	}
	return flux;
}

void Observable::write(IOFiles& w) const {
	w<<links_<<val_;
}

void Observable::print() const {
	if(val_.size()){
		for(unsigned int i(0);i<links_.row();i++){
			std::cout<<links_(i,0)<<" "<<links_(i,1)<<" "<<links_(i,2)<<" "<<val_[links_(i,2)]<<std::endl;
		}
	} else { std::cout<<links_<<std::endl; }
}
/*}*/
