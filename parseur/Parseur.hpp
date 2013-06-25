#ifndef DEF_PARSEUR
#define DEF_PARSEUR

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

class Parseur{
	public:
		/*!Constructor
		 * \param argc simply the argc argument used in int main(int argc, char* argv[])
		 * \param argv[] simply the argc[] argument used in int main(int argc, char* argv[]) */
		Parseur(unsigned int argc, char* argv[]);
		/*! Destroys the static arrays and gives a warning when an argument has not been set */
		~Parseur();

		/*! sets val to the value that corresponds to pattern in argv[]*/
		template<typename Type>
			void set(std::string pattern, Type &input);
		/*! returns the value that corresponds to pattern in argv[]*/
		template<typename Type>
			Type get(std::string pattern);

		bool status() const {return locked; }

	private:
		Parseur();
		Parseur(Parseur const& P);
		Parseur& operator=(Parseur const& P);

		template<typename Type>
			bool search(std::string pattern, Type &input);

		std::vector<std::string> var; //!< vector that contains the options (variables) to set
		std::vector<std::string> val; //!< vector that contains the values of an option
		std::vector<bool> unused; //!< vector that contains false if the corresponding value is used
		bool locked; //!< stores the state of the program, if true a wrong number of agrument was given to the program 
};

template<typename Type>
bool Parseur::search(std::string pattern, Type &input){
	if(!locked){
		unsigned int i(0);
		bool found(false);
		while(!found && i<var.size()){
			if(var[i]==pattern){found=true;}
			else { i++; }
		}
		if(found){
			std::stringstream ss(val[i]);
			ss>>input;
			unused[i] = false;
			return true;

		} else {
			return false;
		}
	} else{
		std::cerr<<"Parseur : the parseur is locked"<<std::endl;
		return true;
	}
}

template<typename Type>
void Parseur::set(std::string pattern, Type &input){
	if(!search(pattern,input)) { std::cerr<<"Parseur : -"<<pattern<<" wasn't found thus its value is unchanged : "<< input <<std::endl; }
}

template<typename Type>
Type Parseur::get(std::string pattern){
	Type input;
	if(!search(pattern, input)) { 
		locked = true; 
		std::cerr<<"Parseur : -"<<pattern<<" wasn't found, the class is locked"<<std::endl; 
	}
	return input;
}
#endif
