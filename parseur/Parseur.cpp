#include"Parseur.hpp"

Parseur::Parseur(unsigned int argc, char* argv[]):
	var(new std::string[argc-1]),
	unused_var(new bool[(argc-1)/2]),
	argc(argc-1),
	locked(false)
{
	if(this->argc % 2 == 0){
		for(unsigned int i(1);i<argc;i++){
			var[i-1] = argv[i];
		}
		for(unsigned int i(1);i<argc/2;i++){
			unused_var[i] = true;
		}
	} else {
		if(this->argc == 1){
			var[0] = argv[1];
			argc=2;
		} else {
			locked = true;
			std::cerr<<"Parseur : not enough arguments"<<std::endl;
		}
	}
}

Parseur::~Parseur(){
	for(unsigned int i(0);i<argc/2;i++){
		if(unused_var[i]){
			std::cerr<<"Parseur : variable "<<var[2*i]<<" was given as input but not used"<<std::endl;
		}
	}
	delete[] var;
	delete[] unused_var;
}
