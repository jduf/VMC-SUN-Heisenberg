#ifndef DEF_TRIANGLEJASTROW
#define DEF_TRIANGLEJASTROW

#include "Triangle.hpp"

class TriangleJastrow: public Triangle<double>{
	public:
		TriangleJastrow(unsigned int N, unsigned int n, unsigned int m);
		~TriangleJastrow();

		void create(double x);
		void check();
		void study();
		void save(Write& w) const;

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

