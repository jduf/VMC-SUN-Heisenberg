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
		/*! returns the value that corresponds to pattern in argv[]*/
		template<typename T>
			T get(std::string pattern);
		/*! work only if argc=2 (one argument passed to the main)*/
		template<typename T>
			void set(T &val);

		unsigned int n_args() const {return argc/2; }

	private:
		Parseur();
		Parseur(Parseur const& P);
		Parseur& operator=(Parseur const& P);

		std::string *var; //!< static array that contains all the char* argv[] strings
		bool *unused_var; //!< static array that lists the variable that where set using void set(std::string pattern, T &val)
		unsigned int argc;//!< size of the arrays (number of arguments given to the main program)
		bool locked; //!< stores the state of the program, if true a wrong number of agrument was given to the program 
};

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

template<typename T>
T Parseur::get(std::string pattern){
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
			T val;
			ss>>val;
			return val;
		} else {
			std::cerr<<"Parseur : "<<pattern<<" wasn't found thus the return value is 0"<<std::endl;
			return 0;
		}
	} else{
		std::cerr<<"Parseur : the parseur is locked because a wrong number of arguments was given"<<std::endl;
		std::cerr<<"Parseur : thus the return value is 0"<<std::endl;
		return 0;
	}
}

template<typename T>
void Parseur::set(T &val){
	if(!locked){
		val = var[0];
	} else{
		std::cerr<<"Parseur : the parseur is locked because a wrong number of arguments was given"<<std::endl;
	}
}
#endif
