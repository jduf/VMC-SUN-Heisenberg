#ifndef DEF_PARSEUR
#define DEF_PARSEUR

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

class Parseur{
	public:
		/*!Constructor
		 * \param argc simply the argc argument used in int main(int argc,
		 * char* argv[])
		 * \param argv[] simply the argc[] argument used in int main(int argc,
		 * char* argv[]) 
		 */
		Parseur(unsigned int argc, char* argv[]);
		Parseur(unsigned int argc, char* argv[], unsigned int N_default);
		/*! Destroys the static arrays and gives a warning when an argument has not been set */
		~Parseur();

		/*! sets val to the value that corresponds to pattern in argv[]*/
		template<typename Type>
			void get(std::string pattern, Type& input);
		/*! returns the value that corresponds to pattern in argv[]*/
		template<typename Type>
			Type get(std::string pattern);

		/*! sets val to the value that corresponds to val[i]*/
		template<typename Type>
			void get(unsigned int i, Type& input);
		/*! returns the value that corresponds to val[i]*/
		template<typename Type>
			Type get(unsigned int i);

		/*! returns true if variable pattern exists*/
		bool check(std::string pattern) const;
		bool status() const {return locked; }

		unsigned int size() const { return val.size();}

		std::string name(unsigned int i) const { return var[i];}

	private:
		Parseur();
		Parseur(Parseur const& P);
		Parseur& operator=(Parseur const& P);

		template<typename Type>
			bool search(std::string pattern, Type &input);

		void init(unsigned int N, char* argv[]);

		std::vector<std::string> var; //!< vector that contains the options (variables) to set
		std::vector<std::string> val; //!< vector that contains the values of an option
		std::vector<bool> used; //!< vector that contains true if the corresponding value is used
		bool locked; //!< stores the state of the program, if true a wrong number of agrument was given to the program 
};

template<typename Type>
bool Parseur::search(std::string pattern, Type &input){
	if(!locked){
		unsigned int i(0);
		while(i<var.size()){
			if(var[i]==pattern){
				std::stringstream ss(val[i]);
				ss>>input;
				used[i] = true;
				return true;
			} else { i++; }
		}
		return false;
	} else {
		std::cerr<<"Parseur : the parseur is locked"<<std::endl;
		return false;
	}
}

template<typename Type>
void Parseur::get(std::string pattern, Type &input){
	if(!search(pattern,input)){ std::cerr<<"Parseur : -"<<pattern<<" wasn't found thus its value is unchanged : "<< input <<std::endl; }
}

template<typename Type>
Type Parseur::get(std::string pattern){
	Type input;
	if(!search(pattern, input)){ 
		locked = true; 
		std::cerr<<"Parseur : -"<<pattern<<" wasn't found, the class is locked"<<std::endl; 
	}
	return input;
}

template<typename Type>
void Parseur::get(unsigned int i, Type &input){
	if(i<val.size()){ input = val[i]; }
	else{std::cerr<<"Parseur : "<<i<<" is bigger than the number of argument given, thus the value is unchanged : "<< input <<std::endl; }
}

template<typename Type>
Type Parseur::get(unsigned int i){
	if(i<val.size()){
		return val[i];
	} else {
		locked = true; 
		std::cerr<<"Parseur : "<<i<<" is bigger than the number of argument given, thus the class is locked"<<std::endl;
		return 0;
	}
}
#endif
