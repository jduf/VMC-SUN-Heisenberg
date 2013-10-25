#ifndef DEF_TRIANGLEMU
#define DEF_TRIANGLEMU

#include "Triangle.hpp"

class TriangleMu: public Triangle<double>{
	public:
		TriangleMu(Parseur& P);
		~TriangleMu();

	protected:
		double mu_;

		void compute_T(unsigned int alpha);
		void save();

		void compute_P();
		void lattice();
		void band_structure();
};
#endif

