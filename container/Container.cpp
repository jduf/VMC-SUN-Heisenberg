#include "Container.hpp"

//Container::Container(bool copyinstring):
	//copyinstring_(copyinstring)
//{}
Container::Container(){}

//Container::Container(Container const& c):
	//copyinstring_(c.copyinstring_),
	//stringvalue_(c.stringvalue_)
//{ 
	//for(unsigned int i(0);i<c.data_.size();i++){
		//data_.push_back((c.data_[i])->clone());
	//}
//}

Container::Container(Container const& c)
	//copyinstring_(c.copyinstring_),
	//stringvalue_(c.stringvalue_)
{ 
	for(unsigned int i(0);i<c.data_.size();i++){
		data_.push_back((c.data_[i])->clone());
	}
}

Container::~Container(){
	for(unsigned int i(0);i<data_.size();i++){
		delete data_[i];
		data_[i] = NULL;
	}
}

//std::string Container::value(unsigned int const& i) const{
	//if(copyinstring_){ return stringvalue_[i]; }
	//else {return ""; }
//}

bool Container::check(std::string pattern) const {
	for(unsigned int i(0);i<data_.size();i++){if(data_[i]->get_name()==pattern){return true;}}
	return false;
}

//std::ostream& operator<<(std::ostream& flux, Container const& c){
	//for(unsigned int i(0);i<c.size();i++){
		//flux<<c.get_as_string(i)<<" ";
	//}
	//return flux;
//}
