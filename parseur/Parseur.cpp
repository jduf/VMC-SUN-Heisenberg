#include "Parseur.hpp"

Parseur::Parseur(unsigned int const& argc, char* argv[]):
	locked_(false)
{
	std::string name;
	std::string val;
	for(unsigned int i(1);i<argc;i++){
		name = argv[i];
		type_.push_back(0);
		if(name[0] == '-'){
			used_.push_back(false);
			/*check if the type is specified*/
			if(name.find(":") != std::string::npos){
				if(i+1<argc){
					val = argv[i+1];
					if(name[1] == 's'){ set(name.substr(3),val); }
					else {
						if(val.find(":") != std::string::npos){ type_.back() = 1; }
						if(val.find(",") != std::string::npos){ type_.back() = 2; }
						switch(type_.back()){
							case 0:
								{
									switch(name[1]){
										case 'i':{ set(name.substr(3),my::string2type<int>(val)); } break;
										case 'u':{ set(name.substr(3),my::string2type<unsigned int>(val)); } break;
										case 'd':{ set(name.substr(3),my::string2type<double>(val)); } break;
										case 'b':{ set(name.substr(3),my::string2type<bool>(val)); } break;
										default: { lock(name); }
									}
								}break;
							case 1:
								{
									switch(name[1]){
										case 'i':{ set_vector_from_range<int>(name,val); } break;
										case 'u':{ set_vector_from_range<unsigned int>(name,val); } break;
										case 'd':{ set_vector_from_range<double>(name,val); } break;
										default: { lock(name); }
									}
								}break;
							case 2:
								{
									switch(name[1]){
										case 'i':{ set_vector_from_list<int>(name,val); } break;
										case 'u':{ set_vector_from_list<unsigned int>(name,val); } break;
										case 'd':{ set_vector_from_list<double>(name,val); } break;
										default: { lock(name); }
									}
								}break;
						}
					}
					i++;
				} else {
					val = "NULL";
					set(name.substr(3),val);
				}
			} else {
				/*!handles option that are linux like (e.g. as in 'ls -l' or
				 * 'head -n 5')*/
				if(i+1<argc){
					/*!if there are other argument in the list, check the
					 * nature of the next one*/
					std::string tmp(argv[i+1]);
					if(tmp[0] == '-'){
						/*!if the next argument has a '-'*/
						if(tmp.find(":") != std::string::npos){
							/*!if the next argument has also a ':' the current
							 * argument can be accepted as a linux like
							 * argument*/
							set(name.substr(1),name.substr(1));
						} else {
							/*!this is problematic because there is no way to
							 * know if it is a value related to the current
							 * argument or a new argument*/
							std::cerr<<__PRETTY_FUNCTION__<<" : problematic argument (impossible to know if '"
								<<tmp<<"' is an option related to '"<<name<<"' or a new argument)"<<std::endl;
							locked_ = true;
						}
					} else {
						/*!in that case the next argument is considered to be
						 * an option related to the current argument*/
						set(name.substr(1),tmp);
						i++;
					}
				} else {
					val = "NULL";
					set(name.substr(1),val);
				}
			}
		} else {
			set(my::tostring(i),name);
			used_.push_back(true);
		}
	}
}

Parseur::~Parseur(){
	for(unsigned int i(0);i<data_.size();i++){
		if(!used_[i]){ std::cerr<<__PRETTY_FUNCTION__<<" : variable "<<data_[i]->name_<<" was given as input but not used"<<std::endl; }
	}
	if(locked_){ std::cerr<<__PRETTY_FUNCTION__<<" : the parseur was locked"<<std::endl; }
}

bool Parseur::find(std::string const& pattern, unsigned int& i, bool lock_iffail) const {
	if(!locked_){
		if(Container::find(pattern,i)){
			used_[i] = true;
			return true;
		} else {
			locked_ = lock_iffail;
			if(locked_){ std::cerr<<__PRETTY_FUNCTION__<<" : can't find "<<pattern<<std::endl; }
			return false;
		}
	} else { return false; }
}

void Parseur::lock(std::string const& arg){
	locked_ = true;
	std::cerr<<__PRETTY_FUNCTION__<<" : wrong argument '"<<arg<<"' : should be '-[iudsb]:name' : "<<std::endl;
}
