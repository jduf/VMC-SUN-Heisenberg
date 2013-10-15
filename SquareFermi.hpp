#ifndef DEF_SQUAREFERMI
#define DEF_SQUAREFERMI

#include "Square.hpp"

class SquareFermi: public Square<double>{
	public:
		SquareFermi(Parseur& P);
		~SquareFermi();

	protected:
		/*! \image html fermi-sea-schema.png*/
		void compute_T();
		void compute_P();
		void compute_band_structure();

		void save(std::string filename);
};
#endif


