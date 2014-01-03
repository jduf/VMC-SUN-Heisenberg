#ifndef DEF_TRIANGLEJASTROW
#define DEF_TRIANGLEJASTROW

#include "Triangle.hpp"
#include "Read.hpp"

class TriangleJastrow: public Triangle<double>{
	public:
		TriangleJastrow(Parseur& P);
		~TriangleJastrow();

		Vector<unsigned int> get_neighbourg(unsigned int i);

	protected:
		Vector<double> nu_;
		Matrix<unsigned int> nn_;
		Vector<unsigned int> sl_;
		Matrix<std::complex<double> > omega_;

		void compute_nn();
		void save();
		void compute_sublattice();
		void compute_omega();
		void lattice();
};
#endif

