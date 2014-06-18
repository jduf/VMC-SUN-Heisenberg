#ifndef DEF_KAGOMEVBC
#define DEF_KAGOMEVBC

#include "Kagome.hpp"
#include "BandStructure.hpp"
#include "PSTricks.hpp"

class KagomeVBC: public Kagome<std::complex<double> >{
	public:
		KagomeVBC(unsigned int const& N, unsigned int const& n, unsigned int const& m, int const& bc, Vector<unsigned int> const& ref);
		~KagomeVBC();

		void create(double const& x, unsigned int const& type);
		void check();
		std::string get_filename() const { return filename_; }

	protected:
		void compute_T();
		void compute_P(Matrix<std::complex<double> >& Px, Matrix<std::complex<double> >& Py);
		void lattice();
};
#endif
