#include"Parseur.hpp"

Parseur::Parseur(unsigned int argc, char* argv[]):
	locked_(false)
{
	if(argv[argc-1][0] == '-'){
		std::cerr<<"warning : Parseur : if the last argument has no var and starts with a -, the prog will crash"<<std::endl;
	}
	unsigned int j(0);
	for(unsigned int i(0);i<argc-1;i++){
		if(argv[i+1][0]=='-'){
			std::string tmp(argv[i+1]);
			var_.push_back(tmp.substr(1));
			val_.push_back(argv[i+2]);
			i++;
		} else {
			std::stringstream ss;
			ss<<j;
			var_.push_back(ss.str());
			val_.push_back(argv[i+1]);
			j++;
		}
		used_.push_back(false);
	}
}

Parseur::~Parseur(){
	for(unsigned int i(0);i<var_.size();i++){
		if(!used_[i]){ std::cerr<<"Parseur : variable "<<var_[i]<<" was given as input but not used"<<std::endl;}
	}
}

bool Parseur::search(std::string const& pattern, unsigned int& i){
	if(!locked_){
		i=0;
		while(i<var_.size()){
			if(var_[i]==pattern){
				return true;
			} else { i++; }
		}
		return false;
	} else {
		std::cerr<<"Parseur : the parseur is locked"<<std::endl;
		return false;
	}
}

std::vector<std::string> &string_split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

std::vector<std::string> string_split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	string_split(s, delim, elems);
	return elems;
}

bool Parseur::is_vector(std::string const& pattern){
	unsigned int i(0);
	if(search(pattern,i)){
		if(val_[i].find(':') != std::string::npos){ return true; } 
		else { return false; } 
	}  else {
		locked_ = true; 
		std::cerr<<"Parseur : -"<<pattern<<" wasn't found, the class is locked"<<std::endl; 
		return false;
	}
}
