#ifndef DEF_SQUAREFERMI
#define DEF_SQUAREFERMI

#include "Square.hpp"
#include "Gnuplot.hpp"

class SquareFermi: public Square<double>{
	public:
		SquareFermi(Parseur& P);
		~SquareFermi();
		void compute_P();

	protected:
		void compute_EVec();
		void compute_spectrum();
		void save(std::string filename);
};
#endif


