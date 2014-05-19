#ifndef DEF_SQUAREJASTROW
#define DEF_SQUAREJASTROW

#include "Square.hpp"

class SquareJastrow: public Square<double>{
	public:
		SquareJastrow(unsigned int N, unsigned int n, unsigned int m);
		~SquareJastrow();

		void create(double x);
		void check();
		void save(IOFiles& w) const;

	protected:
		Matrix<unsigned int> nn_;
		Matrix<unsigned int> cc_;
		Vector<unsigned int> sl_;
		Matrix<std::complex<double> > omega_;

		void compute_nn();
		void compute_sublattice();
		void compute_omega_cc();
		void lattice(Matrix<unsigned int> const& lat);
};
#endif

