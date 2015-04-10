#include"Parseur.hpp"

Parseur::Parseur(unsigned int const& argc, char* argv[]):
	locked_(false)
{
	std::string tmp;
	for(unsigned int i(1);i<argc;i+=2){
		tmp = argv[i];
		if(!locked_ && tmp[0] == '-'){
			if(tmp.size() > 3){
				switch(tmp[2]){
					case ':':
						{
							switch(tmp[1]){
								case 'i':{ set(tmp.substr(3),string2type<int>(argv[i+1])); } break;
								case 'u':{ set(tmp.substr(3),string2type<unsigned int>(argv[i+1])); } break;
								case 'd':{ set(tmp.substr(3),string2type<double>(argv[i+1])); } break;
								case 's':{ set(tmp.substr(3),std::string(argv[i+1])); } break;
								default: { std::cerr<<"Parseur::Parseur(unsigned int argc, char* argv[]) : undefined type"<<std::endl; }
							}
						} break;
					case '{':
						{
							if(tmp[1] == 'd' && tmp[tmp.size()-1] == '}'){
								double min(string2type<double>(argv[i+1]));
								double max(string2type<double>(argv[i+2]));
								double dx (string2type<double>(argv[i+3]));
								std::vector<double> v((max-min)/dx+1);
								for(unsigned int i(0);i<v.size();i++){ v[i] = min+i*dx; }
								set(tmp.substr(3,tmp.size()-4),v);
								is_vec.push_back(size()-1);
								i+=2;
							} else { locked_ = true; }
						} break;
					default: { locked_ = true; }
				}
				used_.push_back(false);
			} else {
				set(tmp.substr(1),std::string(argv[i+1])); 
				used_.push_back(true);
			}
		} else {
			i--;
			set(my::tostring(i),std::string(argv[i+1])); 
			used_.push_back(true);
		}
	}
	if(locked_){
		std::cerr<<"Parseur::Parseur(unsigned int argc, char* argv[]) : wrong argument, should be '-[iuds]:name' or '-[uids]{name}'"<<std::endl; 
	}
}

Parseur::~Parseur(){
	for(unsigned int i(0);i<data_.size();i++){
		if(!used_[i]){ std::cerr<<"Parseur : variable "<<data_[i]->get_name()<<" was given as input but not used"<<std::endl;}
	}
	if(locked_){ std::cerr<<"Parseur::~Parseur() : the parseur was locked"<<std::endl; }
}

bool Parseur::find(std::string const& pattern, unsigned int& i, bool iffail){
	if(!locked_){
		i=0;
		if(Container::find(pattern,i)){
			used_[i] = true;
			return true;
		} else {
			locked_ = iffail;
			return false;
		}
	} else {
		return false;
	}
}

bool Parseur::is_vector(std::string const& pattern){
	unsigned int i(0);
	if(!locked_ && find(pattern,i)){
		unsigned int j(0);
		while(j<is_vec.size()){ 
			if(i == is_vec[j] ) {return true;}
			j++;
		}
		return false;
	}
	locked_ = true; 
	std::cerr<<"bool Parseur::is_vector(std::string const& pattern) : -"<<pattern<<" wasn't found, the class is locked"<<std::endl; 
	return false;
}
