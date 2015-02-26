#ifndef DEF_PARSEUR
#define DEF_PARSEUR

#include "Container.hpp"
#include "Vector.hpp"

class Parseur: public Container {
	public:
		/*{Description*/
		/*!Constructor
		 * \param argc=argc argument used in main(int argc, char* argv[])
		 * \param argv[]=argv[] argument used in main(int argc, char* argv[]) 
		 */
		/*}*/
		Parseur(unsigned int const& argc, char* argv[]);
		/*!Destructor that warns when an argument has not been set*/
		~Parseur();

		/*!Returns true if there is a var_[i]==patern and sets i*/
		bool find(std::string const& pattern, unsigned int& i, bool iffail=true);
		/*!Returns true if val_[i] for var_[i]==patern is a vector*/
		bool is_vector(std::string const& pattern);
		/*!Returns locked_*/
		bool status() const { return locked_; }

	private:
		/*!Forbids default*/
		Parseur();
		/*!Forbids copy*/
		Parseur(Parseur const& P);
		/*!Forbids assignment*/
		Parseur& operator=(Parseur const& P);

		/*!Uses std::stringstream to set input=val_[i]*/
		template<typename Type>
			Type string2type(std::string const& s);
		
		std::vector<bool> used_;
		std::vector<unsigned int> is_vec;
		bool locked_;
};

std::vector<std::string> &string_split(const std::string &s, char delim, std::vector<std::string> &elems);

std::vector<std::string> string_split(const std::string &s, char delim);

template<typename Type>
Type Parseur::string2type(std::string const& s){
	Type out;
	std::stringstream ss(s);
	ss>>out;
	return out;
}
#endif
