#include "Observable.hpp"

Observable::Observable(unsigned int const& n, unsigned int const& B, unsigned int const& b, bool const& conv):
	links_(n,2)
{
	val_.set(n,B,b,conv);
}

Observable::Observable(IOFiles& r):
	links_(r),
	val_(r)
{}

IOFiles& operator<<(IOFiles& w, Observable const& obs){
	if(w.is_binary()){
		(void)(obs);
	}
	return w;
}

IOFiles& operator>>(IOFiles& r, Observable& obs){
	if(r.is_binary()){
		(void)(obs);
	} 
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
