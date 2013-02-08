#ifndef DEF_PARSEUR
#define DEF_PARSEUR

#include <iostream>
#include <string>
#include <sstream>

class Parseur{
	public:
		Parseur(unsigned int argc, char* argv[]);
		~Parseur();

		template<typename T>
		void set(std::string pattern, T &val);

	private:
		std::string *var;
		unsigned int argc;
};

Parseur::Parseur(unsigned int argc, char* argv[]):
	var(new std::string[argc-1]),
	argc(argc-1)
{
	if(this->argc % 2 == 0){
		for(unsigned int i(1);i<argc;i++){
			var[i-1] = argv[i];
		}
	} else {
		std::cerr<<"Parseur : not enough arguments"<<std::endl;
	}
}

Parseur::~Parseur(){
	delete[] var;
}

template<typename T>
void Parseur::set(std::string pattern, T &val){
	bool found(false);
	unsigned int i(0);
	while(!found && i<argc){
		if(var[i]==pattern){found = true;}
		else{i += 2;}
	}
	if(found){
		std::stringstream ss(var[i+1]);
		ss>> val;
	} else {
		std::cerr<<"Parseur : "<<pattern<<" wasn't found thus its value is not defined"<<std::endl;
	}
}
#endif
