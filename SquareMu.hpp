#ifndef DEF_SQUAREMU
#define DEF_SQUAREMU

#include "Square.hpp"

class SquareMu: public Square<double>{
	public:
		SquareMu(Parseur& P);
		~SquareMu();

	protected:
		double mu;

		void compute_EVec(unsigned int spin);
		void save(std::string filename);
		void show(Matrix<double> const& T,unsigned int spin);
};
#endif

