#ifndef DEF_SQUAREJASTROW
#define DEF_SQUAREJASTROW

#include "Square.hpp"

class SquareJastrow: public Square<double>{
	public:
		SquareJastrow(Parseur& P);
		~SquareJastrow();

	protected:
		Vector<double> nu_;
		Matrix<unsigned int> nn_;
		Matrix<unsigned int> nnn_;
		Vector<unsigned int> sl_;
		Matrix<std::complex<double> > omega_;

		Vector<unsigned int> get_neighbourg(unsigned int i);
		void compute_nn();
		void compute_nnn();
		void save();
		void compute_sublattice();
		void compute_omega();
};
#endif

