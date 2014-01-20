#ifndef DEF_SQUAREJASTROW
#define DEF_SQUAREJASTROW

#include "Square.hpp"
#include "Container.hpp"

class SquareJastrow: public Square<double>{
	public:
		SquareJastrow(Parseur& P);
		~SquareJastrow();

		void properties(Container& c);
		void save();
		void lattice(Matrix<unsigned int> const& lat);

	protected:
		Matrix<unsigned int> nn_;
		Matrix<unsigned int> cc_;
		Vector<unsigned int> sl_;
		Matrix<std::complex<double> > omega_;

		void compute_nn();
		void compute_sublattice();
		void compute_omega_cc();
};
#endif

