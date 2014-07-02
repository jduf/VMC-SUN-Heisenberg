#ifndef DEF_PARSEUR
#define DEF_PARSEUR

#include <vector>
#include "Vector.hpp"

class Parseur{
	public:
		/*{Description*/
		/*!Constructor
		 * \param argc=argc argument used in main(int argc, char* argv[])
		 * \param argv[]=argv[] argument used in main(int argc, char* argv[]) 
		 */
		/*}*/
		Parseur(unsigned int argc, char* argv[]);
		/*!Destructor that warns when an argument has not been set */
		~Parseur();

		/*!Sets input to the val_[i] such that var_[i]==pattern*/
		template<typename Type>
			void get(std::string const& pattern, Type& input);
		/*!Returns the val_[i] such that var_[i]==pattern*/
		template<typename Type>
			Type get(std::string const& pattern);

		/*!Sets input to val_[i]*/
		template<typename Type>
			void get(unsigned int i, Type& input);
		/*!Returns val_[i]*/
		template<typename Type>
			Type get(unsigned int i);

		/*!Returns locked_*/
		bool status() const {return locked_; }
		/*!Returns true if there is a var_[i]==patern and sets i*/
		bool search(std::string const& pattern, unsigned int& i);
		/*!Returns true if val_[i] for var_[i]==patern is a vector*/
		bool is_vector(std::string const& pattern);

		/*!Returns the number of values stored in the class*/
		unsigned int size() const { return val_.size();}
		/*!Returns the name of the variable i*/
		std::string name(unsigned int i) const { return var_[i];}

	private:
		/*!Forbids default*/
		Parseur();
		/*!Forbids copy*/
		Parseur(Parseur const& P);
		/*!Forbids assignment*/
		Parseur& operator=(Parseur const& P);

		/*!Uses std::stringstream to set input=val_[i]*/
		template<typename Type>
			void set(unsigned int i, Type &input);
		/*!Uses std::stringstream to set input=val_[i] for a Vector*/
		template<typename Type>
			void set(unsigned int i, Vector<Type> &input);

		std::vector<std::string> var_; //!< vector that contains the options (variables) to set
		std::vector<std::string> val_; //!< vector that contains the values of an option
		std::vector<bool> used_; //!< used_[i]==true if val_[i] has been used
		bool locked_; //!< true if a wrong number of agrument was given to the program 
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
		locked_ = true; 
		std::cerr<<"Parseur : -"<<pattern<<" wasn't found, the class is locked"<<std::endl; 
		return Type();
	}
}

template<typename Type>
void Parseur::get(unsigned int i, Type &input){
	if(i<val_.size()){set(i,input);}
	else{ std::cerr<<"Parseur : the "<<i<<"th entry wasn't found thus its value is unchanged : "<< input <<std::endl; }
}

template<typename Type>
Type Parseur::get(unsigned int i){
	if(i<val_.size()){
		Type input;
		set(i,input);
		return input;
	} else {
		locked_ = true; 
		std::cerr<<"Parseur : the "<<i<<"th entry wasn't found, the class is locked"<<std::endl; 
		return Type();
	}
}

template<typename Type>
void Parseur::set(unsigned int i, Type &input){
	std::stringstream ss(val_[i]);
	ss>>input;
	used_[i] = true;
}

template<typename Type>
void Parseur::set(unsigned int i, Vector<Type> &input){
	std::vector<std::string> v(string_split(val_[i],':'));
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
	for(unsigned int j(0);j<n;j++){ input(j) = a+j*d; }
	used_[i] = true;
}
#endif
