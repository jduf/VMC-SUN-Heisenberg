#include"Parseur.hpp"

Parseur::Parseur(unsigned int argc, char* argv[]):
	locked(false)
{
	if(argv[argc-1][0] == '-'){
		std::cerr<<"warning : Parseur : if the last argument has no var and starts with a -, the prog will crash"<<std::endl;
	}
	init(argc-1, argv);
}

Parseur::Parseur(unsigned int argc, char* argv[], unsigned int N_default):
	locked(false)
{
	if( 2*N_default < argc ) { 
		if(argv[argc-1][0] == '-'){
			std::cerr<<"warning : Parseur : if the last argument has no var and starts with a -, the prog will crash"<<std::endl;
		}
		init(2*N_default,argv);

		for(unsigned int i(2*N_default+1);i<argc;i++){
			val.push_back(argv[i]);
		}
	} else {
		locked=true;
		std::cerr<<"warning : Parseur : you should give "<< N_default<<" arguments"<<std::endl;
	}
}

Parseur::~Parseur(){
	for(unsigned int i(0);i<var.size();i++){
		if(!used[i]){ std::cerr<<"Parseur : variable "<<var[i]<<" was given as input but not used"<<std::endl;}
	}
}

void Parseur::init(unsigned int N, char* argv[]){
	unsigned int j(0);
	for(unsigned int i(0);i<N;i++){
		if(argv[i+1][0]=='-'){
			std::string tmp(argv[i+1]);
			var.push_back(tmp.substr(1));
			val.push_back(argv[i+2]);
			i++;
		} else {
			std::stringstream ss;
			ss<<j;
			var.push_back(ss.str());
			val.push_back(argv[i+1]);
			j++;
		}
		used.push_back(false);
	}
}

bool Parseur::search(std::string const& pattern, unsigned int& i){
	if(!locked){
		i=0;
		while(i<var.size()){
			if(var[i]==pattern){
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
		if(val[i].find(':') != std::string::npos){ return true; } 
		else { return false; } 
	}  else {
		locked = true; 
		std::cerr<<"Parseur : -"<<pattern<<" wasn't found, the class is locked"<<std::endl; 
		return false;
	}
}
