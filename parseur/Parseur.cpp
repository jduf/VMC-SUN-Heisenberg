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
								set(tmp.substr(3,tmp.size()-4), 
										Vector<double>(string2type<double>(argv[i+1]),
											string2type<double>(argv[i+2]),
											string2type<double>(argv[i+3])));
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
			set(tostring(i),std::string(argv[i+1])); 
			used_.push_back(true);
		}
	}
	if(locked_){
		std::cerr<<"Parseur::Parseur(unsigned int argc, char* argv[]) : wrong argument, should be '-[iuds]:name'"<<std::endl; 
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
			std::cout<<pattern<<std::endl;
			return true;
		} else {
			std::cout<<"not"<<pattern<<std::endl;
			locked_ = iffail;
			return false;
		}
	} else {
		return false;
	}
}

std::vector<std::string> &string_split(const std::string &s, char delim, std::vector<std::string> &elems){
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) { elems.push_back(item); }
	return elems;
}

std::vector<std::string> string_split(const std::string &s, char delim){
	std::vector<std::string> elems;
	string_split(s, delim, elems);
	return elems;
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
	std::cerr<<"Parseur : -"<<pattern<<" wasn't found, the class is locked"<<std::endl; 
	return false;
}
