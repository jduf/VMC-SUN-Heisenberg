#ifndef DEF_TRIANGLEJASTROW
#define DEF_TRIANGLEJASTROW

#include "Triangle.hpp"
#include "Container.hpp"

class TriangleJastrow: public Triangle<double>{
	public:
		TriangleJastrow(Parseur& P);
		~TriangleJastrow();

		void properties(Container& c);
		void save();
		void lattice(Matrix<unsigned int> const& lat);

	protected:
		Matrix<unsigned int> nn_;
		Matrix<unsigned int> cc_;
		Vector<unsigned int> sl_;
		Matrix<std::complex<double> > omega_;

		Vector<unsigned int> get_neighbourg(unsigned int i);
		void compute_nn();
		void compute_sublattice();
		void compute_omega_cc();
};
#endif

