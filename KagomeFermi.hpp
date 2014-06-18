#ifndef DEF_KAGOMEFERMI
#define DEF_KAGOMEFERMI

#include "Kagome.hpp"
#include "BandStructure.hpp"

class KagomeFermi: public Kagome<double>{
	public:
		KagomeFermi(unsigned int const& N, unsigned int const& n, unsigned int const& m, int const& bc, Vector<unsigned int> const& ref);
		~KagomeFermi();

		void create(double const& x, unsigned int const& type);
		void check();
		std::string get_filename() const { return filename_; }

	protected:
		void compute_T();
		void compute_P(Matrix<double>& Px, Matrix<double>& Py);
		void lattice();
};
#endif
