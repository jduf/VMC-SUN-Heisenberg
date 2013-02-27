#ifndef DEF_PARSEUR
#define DEF_PARSEUR

#include <iostream>
#include <string>
#include <sstream>

class Parseur{
	public:
		/*!Constructor
		 * \param argc simply the argc argument used in int main(int argc, char* argv[])
		 * \param argv[] simply the argc[] argument used in int main(int argc, char* argv[]) */
		Parseur(unsigned int argc, char* argv[]);
		/*! Destroys the static arrays and gives a warning when an argument has not been set */
		~Parseur();

		/*! sets val to the value that corresponds to pattern in argv[]*/
		template<typename T>
		void set(std::string pattern, T &val);

	private:
		Parseur();
		Parseur(Parseur const& P);
		Parseur& operator=(Parseur const& P);

		std::string *var; //!< static array that contains all the char* argv[] strings
		bool *unused_var; //!< static array that lists the variable that where set using void set(std::string pattern, T &val)
		unsigned int argc;//!< size of the arrays (number of arguments given to the main program 
		bool locked; //!< stores the state of the program, if true a wrong number of agrument was given to the program 
};

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
		locked = true;
		std::cerr<<"Parseur : not enough arguments"<<std::endl;
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

template<typename T>
void Parseur::set(std::string pattern, T &val){
	if(!locked){
		bool found(false);
		unsigned int i(0);
		while(!found && i<argc){
			if(var[i]==pattern){found = true;}
			else{i += 2;}
		}
		if(found){
			unused_var[i/2] = false;
			std::stringstream ss(var[i+1]);
			ss>>val;
		} else {
			std::cerr<<"Parseur : "<<pattern<<" wasn't found thus its value is still "<<val<<std::endl;
		}
	} else{
		std::cerr<<"Parseur : the parseur is locked because a wrong number of arguments was given"<<std::endl;
	}
}
#endif
