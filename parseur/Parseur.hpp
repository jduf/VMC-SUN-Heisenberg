#ifndef DEF_PARSEUR
#define DEF_PARSEUR

#include "Container.hpp"

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
		/*{Forbidden*/
		Parseur() = delete;
		Parseur(Parseur const&) = delete;
		Parseur(Parseur&&) = delete;
		Parseur& operator=(Parseur) = delete;
		/*}*/

		/*!Returns true if there is a var_[i]==patern and sets i*/
		bool find(std::string const& pattern, unsigned int& i, bool iffail=true);
		/*!Returns true if val_[i] for var_[i]==patern is a vector*/
		bool is_vector(std::string const& pattern);
		/*!Returns locked_*/
		bool status() const { return locked_; }

	private:
		/*!Uses std::stringstream to set input=val_[i]*/
		template<typename Type>
			Type string2type(std::string const& s);
		
		std::vector<bool> used_;
		std::vector<unsigned int> is_vec;
		bool locked_;
};

template<typename Type>
Type Parseur::string2type(std::string const& s){
	Type out;
	std::stringstream ss(s);
	ss>>out;
	return out;
}
#endif
