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
		/*! work only if argc=2 (one argument passed to the main)*/
		//template<typename Type>
			//void set(Type &val);
		/*! returns the value that corresponds to pattern in argv[]*/
		template<typename Type>
			Type get(std::string pattern);

		unsigned int n_args() const {return var.size(); }
		bool status() const {return locked; }

		void print();

	private:
		Parseur();
		Parseur(Parseur const& P);
		Parseur& operator=(Parseur const& P);

		std::vector<std::string> var; //!< vector that contains the options (variables) to set
		std::vector<std::string> val; //!< vector that contains the values of an option
		bool locked; //!< stores the state of the program, if true a wrong number of agrument was given to the program 
};

template<typename Type>
void Parseur::set(std::string pattern, Type &input){
	if(!locked){
		bool found(false);
		unsigned int i(0);
		while(!found && i<var.size()){
			if(var[i]==pattern){found = true;}
			else{i++;}
		}
		if(found){
			std::stringstream ss(val[i]);
			ss>>input;
			val.erase(val.begin()+i);
			var.erase(var.begin()+i);
		} else {
			std::cerr<<"Parseur : -"<<pattern<<" wasn't found"<<std::endl;
			locked=true;
		}
	} else{
		std::cerr<<"Parseur : the parseur is locked because a wrong number of arguments was given"<<std::endl;
	}
}

template<typename Type>
Type Parseur::get(std::string pattern){
	Type input;
	set(pattern,input);
	return input;
}

//template<typename Type>
//void Parseur::set(Type &val){
	//if(!locked){
		//val = var[0];
	//} else{
		//std::cerr<<"Parseur : the parseur is locked because a wrong number of arguments was given"<<std::endl;
	//}
//}
#endif
