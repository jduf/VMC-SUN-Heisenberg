#include "Container.hpp"

Data::~Data(){}

Container::Container(bool copyinstring):
	copyinstring_(copyinstring)
{ }

Container::Container(Container const& c):
	copyinstring_(c.copyinstring_),
	stringvalue_(c.stringvalue_)
{ 
	for(unsigned int i(0);i<c.d_.size();i++){
		d_.push_back((c.d_[i])->clone());
	}
}

Container::~Container(){
	for(unsigned int i(0);i<d_.size();i++){
		delete d_[i];
		d_[i] = NULL;
	}
}

std::string Container::value(unsigned int const& i) const{
	if(copyinstring_){ return stringvalue_[i]; }
	else {return ""; }
}

std::ostream& operator<<(std::ostream& flux, Container const& input){
	for(unsigned int i(0);i<input.size();i++){
		flux<<input.value(i)<<" ";
	}
	return flux;
}

bool Container::check(std::string pattern) const {
	for(unsigned int i(0);i<d_.size();i++){if(d_[i]->get_name()==pattern){return true;}}
	return false;
}
