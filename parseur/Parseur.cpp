#include"Parseur.hpp"

Parseur::Parseur(unsigned int const& argc, char* argv[]):
	locked_(false)
{
	std::string name;
	std::string val;
	unsigned int type;
	for(unsigned int i(1);i<argc;i+=2){
		name = argv[i];
		used_.push_back(false);
		if(name[0] == '-'){
			val = argv[i+1];
			/*check if the type is specified*/
			if(name.find(":") != std::string::npos){
				if(name[1] == 's'){
					set(name.substr(3),std::string(val));
				} else{
					type = 0;
					if(val.find(":") != std::string::npos){ type = 1; }
					if(val.find(",") != std::string::npos){ type = 2; }
					switch(type){
						case 0:
							{
								switch(name[1]){
									case 'i':{ set(name.substr(3),my::string2type<int>(val)); } break;
									case 'u':{ set(name.substr(3),my::string2type<unsigned int>(val)); } break;
									case 'd':{ set(name.substr(3),my::string2type<double>(val)); } break;
									default: { lock(name); }
								}
							}break;
						case 1:
							{
								switch(name[1]){
									case 'i':{set_vector_from_range<int>(name,val);} break;
									case 'u':{set_vector_from_range<unsigned int>(name,val);} break;
									case 'd':{set_vector_from_range<double>(name,val);} break;
									default: { lock(name); }
								}
							}break;
						case 2:
							{
								switch(name[1]){
									case 'i':{set_vector_from_list<int>(name,val);} break;
									case 'u':{set_vector_from_list<unsigned int>(name,val);} break;
									case 'd':{set_vector_from_list<double>(name,val);} break;
									default: { lock(name); }
								}
							}break;
					}
				}
			} else { set(name.substr(1),val); }
		} else {
			set(my::tostring(i-1),std::string(argv[i])); 
			i--;
		}
	}
	if(locked_){
	}
}

Parseur::~Parseur(){
	for(unsigned int i(0);i<data_.size();i++){
		if(!used_[i]){ std::cerr<<"Parseur : variable "<<data_[i]->get_name()<<" was given as input but not used"<<std::endl;}
	}
	if(locked_){ std::cerr<<"Parseur::~Parseur() : the parseur was locked"<<std::endl; }
}

bool Parseur::find(std::string const& pattern, unsigned int& i, bool iffail) const {
	if(!locked_){
		if(Container::find(pattern,i)){
			used_[i] = true;
			return true;
		} else {
			locked_ = iffail;
			return false;
		}
	} else { return false; }
}

void Parseur::lock(std::string const& arg){
	locked_ = true;
	std::cerr<<"Parseur::Parseur(unsigned int argc, char* argv[]) : wrong argument '"<<arg<<"' : should be '-[iuds]:name' or '-[uids].name' : "<<std::endl; 
}

