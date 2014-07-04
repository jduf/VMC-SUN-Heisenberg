#ifndef DEF_HONEYCOMB0PP
#define DEF_HONEYCOMB0PP

#include "Honeycomb.hpp"

class Honeycomb0pp: public Honeycomb<double>{
	public:
		Honeycomb0pp(Vector<unsigned int> const& ref, unsigned int const& N, unsigned int const& m, unsigned int const& n, Vector<unsigned int> const& M,  int const& bc, double td);
		~Honeycomb0pp(){}

		void create();
		void save(IOFiles& w) const;
		void check();

	protected:
		double td_;

		void compute_T();
		void lattice();
		void compute_P(Matrix<double>& Px, Matrix<double>& Py);
};
#endif

