#ifndef DEF_TRIANGLEMU
#define DEF_TRIANGLEMU

#include "Triangle.hpp"

class TriangleMu: public Triangle<double>{
	public:
		TriangleMu(Parseur& P);
		~TriangleMu();

		void save();

	protected:
		double mu_;

		void compute_T(unsigned int alpha);
		void compute_P();
		void lattice();
};
#endif

