#ifndef DEF_SQUAREFERMI
#define DEF_SQUAREFERMI

#include "Square.hpp"

class SquareFermi: public Square<double>{
	public:
		SquareFermi(Parseur& P);
		~SquareFermi();

	protected:
		void compute_T();
		void save();

		void compute_P();
		void lattice();
		void band_structure();
};
#endif
