#ifndef DEF_SQUAREPIFLUX
#define DEF_SQUAREPIFLUX

#include "Square.hpp"

class SquarePiFlux: public Square<std::complex<double> >{
	public:
		SquarePiFlux(unsigned int const& N, unsigned int const& n, unsigned int const& m, int const& bc, Vector<unsigned int> const& ref);
		~SquarePiFlux();

		void create(double const& x, unsigned int const& type);
		void check();
		std::string get_filename() const { return filename_;}

		void treat_one_sim(IOFiles& read, IOFiles& write, RSTFile& rst, std::string const& path, std::string const& filename){}

	protected:
		void compute_T();
};
#endif
