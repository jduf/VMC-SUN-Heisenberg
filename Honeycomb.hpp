#ifndef DEF_HONEYCOMB
#define DEF_HONEYCOMB

#include "CreateSystem.hpp"

class Honeycomb: public CreateSystem<double>{
	public:
		Honeycomb(Parseur& P);
		~Honeycomb();

	protected:
		unsigned int N_row, N_col;

		void compute_T();
};
#endif

