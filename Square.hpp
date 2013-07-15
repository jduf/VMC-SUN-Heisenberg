#ifndef DEF_SQUARE
#define DEF_SQUARE

#include "CreateSystem.hpp"

#include <complex>

class Square: public CreateSystem<std::complex<double> >{
	public:
		Square(Parseur& P);
		~Square();

	protected:
		unsigned int N_row, N_col;

		void compute_T();
		void save(std::string filename);
};

#endif

