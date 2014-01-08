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

	protected:
		Matrix<unsigned int> nn_;
		Vector<unsigned int> sl_;
		Matrix<std::complex<double> > omega_;

		Vector<unsigned int> get_neighbourg(unsigned int i);
		void compute_nn();
		void compute_sublattice();
		void compute_omega();
};
#endif

