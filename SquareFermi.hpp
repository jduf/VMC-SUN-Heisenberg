#ifndef DEF_SQUAREFERMI
#define DEF_SQUAREFERMI

#include "Square.hpp"

class SquareFermi: public Square<double>{
	public:
		SquareFermi(Parseur& P);
		~SquareFermi();

	protected:
		void compute_T();
		void compute_P();
		void save(std::string filename);
};
#endif


