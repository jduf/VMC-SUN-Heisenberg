#ifndef DEF_SQUAREMU
#define DEF_SQUAREMU

#include "Square.hpp"

class SquareMu: public Square<double>{
	public:
		SquareMu(Parseur& P);
		~SquareMu();

	protected:
		double mu_;

		/*! \image html stripe-order-brillouin.png*/
		void compute_T(unsigned int spin);
		void compute_P();
		void compute_band_structure(unsigned int spin);

		void save(std::string filename);
		void show(Matrix<double> const& T,unsigned int spin);
};
#endif

