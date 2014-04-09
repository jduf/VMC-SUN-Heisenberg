#ifndef DEF_PARSEUR
#define DEF_PARSEUR

#include <string>
#include <sstream>
#include <vector>

#include "Vector.hpp"

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
			void get(std::string const& pattern, Type& input);
		/*! returns the value that corresponds to pattern in argv[]*/
		template<typename Type>
			Type get(std::string const& pattern);

		/*! sets val to the value that corresponds to val[i]*/
		template<typename Type>
			void get(unsigned int i, Type& input);
		/*! returns the value that corresponds to val[i]*/
		template<typename Type>
			Type get(unsigned int i);

		bool status() const {return locked; }

		bool search(std::string const& pattern, unsigned int& i);
		bool is_vector(std::string const& pattern);

		unsigned int size() const { return val.size();}

		std::string name(unsigned int i) const { return var[i];}

	private:
		Parseur();
		Parseur(Parseur const& P);
		Parseur& operator=(Parseur const& P);

		template<typename Type>
			void set(unsigned int i, Type &input);
		template<typename Type>
			void set(unsigned int i, Vector<Type> &input);

		void init(unsigned int N, char* argv[]);

		std::vector<std::string> var; //!< vector that contains the options (variables) to set
		std::vector<std::string> val; //!< vector that contains the values of an option
		std::vector<bool> used; //!< vector that contains true if the corresponding value is used
		bool locked; //!< stores the state of the program, if true a wrong number of agrument was given to the program 
};

std::vector<std::string> &string_split(const std::string &s, char delim, std::vector<std::string> &elems);

std::vector<std::string> string_split(const std::string &s, char delim);

template<typename Type>
void Parseur::get(std::string const& pattern, Type &input){
	unsigned int i(0);
	if(search(pattern,i)){set(i,input);}
	else{ std::cerr<<"Parseur : -"<<pattern<<" wasn't found thus its value is unchanged : "<< input <<std::endl; }
}

template<typename Type>
Type Parseur::get(std::string const& pattern){
	unsigned int i(0);
	if(search(pattern,i)){ 
		Type input;
		set(i,input);
		return input;
	} else {
		locked = true; 
		std::cerr<<"Parseur : -"<<pattern<<" wasn't found, the class is locked"<<std::endl; 
		return Type();
	}
}

template<typename Type>
void Parseur::get(unsigned int i, Type &input){
	if(i<val.size()){set(i,input);}
	else{ std::cerr<<"Parseur : the "<<i<<"th entry wasn't found thus its value is unchanged : "<< input <<std::endl; }
}

template<typename Type>
Type Parseur::get(unsigned int i){
	if(i<val.size()){
		Type input;
		set(i,input);
		return input;
	} else {
		locked = true; 
		std::cerr<<"Parseur : the "<<i<<"th entry wasn't found, the class is locked"<<std::endl; 
		return Type();
	}
}

template<typename Type>
void Parseur::set(unsigned int i, Type &input){
	std::stringstream ss(val[i]);
	ss>>input;
	used[i] = true;
}

template<typename Type>
void Parseur::set(unsigned int i, Vector<Type> &input){
	std::vector<std::string> v(string_split(val[i],':'));
	std::stringstream sa(v[0]);
	std::stringstream sb(v[2]);
	std::stringstream sd(v[1]);
	double a;
	double b;
	double d;
	sa>>a;
	sb>>b;
	sd>>d;
	unsigned int n((b-a)/d+1);
	input.set(n);
	for(unsigned int j(0);j<n;j++){
		input(j) = a+j*d;
	}
	used[i] = true;
}
#endif
