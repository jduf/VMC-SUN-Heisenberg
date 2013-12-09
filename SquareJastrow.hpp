#ifndef DEF_SQUAREJASTROW
#define DEF_SQUAREJASTROW

#include "Square.hpp"

class SquareJastrow: public Square<double>{
	public:
		SquareJastrow(Parseur& P);
		~SquareJastrow();

	protected:
		double nu_;

		//void compute_T();
		void save();

		//void compute_P();
		//void band_structure();
		//void lattice();
};
#endif

