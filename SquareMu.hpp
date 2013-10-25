#ifndef DEF_SQUAREMU
#define DEF_SQUAREMU

#include "Square.hpp"

class SquareMu: public Square<double>{
	public:
		SquareMu(Parseur& P);
		~SquareMu();

	protected:
		double mu_;

		void compute_T(unsigned int alpha);
		void save();

		void compute_P();
		void lattice();
		void band_structure();
};
#endif

