#ifndef DEF_CHAINFERMI
#define DEF_CHAINFERMI

#include "Chain.hpp"

class ChainFermi: public Chain<double>{
	public:
		ChainFermi(unsigned int const& N, unsigned int const& n, unsigned int const& m, int const& bc, Vector<unsigned int> const& ref);
		~ChainFermi();

		void create(double const& delta, unsigned int const& type);
		void check();
		std::string get_filename() const { return filename_;}

	private:
		void compute_P(Matrix<double>& P);
		void compute_T();
};
#endif
